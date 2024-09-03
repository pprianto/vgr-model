function make_Variables(
    model::Model,
    sets::NamedTuple,
    params::NamedTuple
)

#=------------------------------------------------------------------------------
------------------------------ VARIABLES ---------------------------------------

Define Variables of the model
current variables:
1. generation and storage capacity
2. generation and storage dispatch
3. power flow
4. voltage level
5. import/export
6. ...

------------------------------------------------------------------------------=#
    ## Sets and parameters

    @unpack NODES, 
            TRANSMISSION_NODES,
            SE3_TRANS_NODES,
            NO1_TRANS_NODES,
            DK1_TRANS_NODES,
            COAST_NODES,
            GBG,
            GEN_TECHS, 
            EL_GEN, 
            HEAT_GEN, 
            H2_GEN, 
            STO_TECHS, 
            EL_STO, 
            HEAT_STO, 
            H2_STO, 
            PERIODS, 
            LINES, 
            NODE_FROM, 
            NODE_TO, 
            CHP,
            FC,
            WIND,
            PV,
            HP,
            BOILER,
            EC,
            FLEX_TH,
            THERMAL_1H,
            THERMAL_2H,
            # THERMAL_8H,
            THERMAL_12H = sets
    
    @unpack SE3_price,
            NO1_price,
            DK1_price,  
            Gentech_data, 
            Stotech_data, 
            Eldemand_data, 
            Heatdemand_data, 
            H2demand_data, 
            Discount_rate,
            Gen_cos_ϕ, 
            Gen_sin_ϕ, 
            Demand_cos_ϕ, 
            Demand_sin_ϕ, 
            Lines_props,
            Vnom = params

    #=------------------------------------------------------------------------------
    MODEL VARIABLES
    ------------------------------------------------------------------------------=#

    #=------------------------------------------------------------------------------
    GENERATION AND STORAGE VARIABLES
    ------------------------------------------------------------------------------=#
    # EQ (2) - (7) #
    # Generation and Storage Capacities (MW)
    # upper and lower bounds for generation capaicity is defined in set_gen_bounds function
    @variables model begin
        existing_generation[i ∈ NODES, x ∈ GEN_TECHS] ≥ 0
    end
    
    @variables model begin
        generation_investment[i ∈ NODES, x ∈ GEN_TECHS]   ≥ 0
        storage_investment[i ∈ NODES, s ∈ STO_TECHS]      ≥ 0
    end
    
    # Active and Reactive power dispatch (MWh or MVArh)
    # reactive generation is assumed only applies for electricity generation technologies
    @variables model begin
        active_generation[i ∈ NODES, x ∈ GEN_TECHS, t ∈ PERIODS]  ≥ 0
        reactive_generation[i ∈ NODES, x ∈ EL_GEN, t ∈ PERIODS]   ≥ 0
    end

    # Storage-related variables (MWh)
    @variables model begin
        storage_charge[i ∈ NODES, s ∈ STO_TECHS, t ∈ PERIODS]     ≥ 0
        storage_discharge[i ∈ NODES, s ∈ STO_TECHS, t ∈ PERIODS]  ≥ 0
        storage_level[i ∈ NODES, s ∈ STO_TECHS, t ∈ PERIODS]      ≥ 0
    end

    #=------------------------------------------------------------------------------
    TWO-VARIABLE APPROACH VARIABLES
    ------------------------------------------------------------------------------=#
    # applies to thermal power plants
    # related to start up / ramping limits
    # based on Lisa's thesis
    @variables model begin
        generation_spin[i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS]      ≥ 0
        generation_on[i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS]        ≥ 0
        gen_startup_cost[i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS]     ≥ 0
        gen_partload_cost[i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS]    ≥ 0
        gen_startup_CO2[i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS]      ≥ 0
        gen_partload_CO2[i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS]     ≥ 0
    end

    #=------------------------------------------------------------------------------
    EXPORT - IMPORT VARIABLES
    EQ (38)
    ------------------------------------------------------------------------------=#
    # EQ (33) - (35) #
    # export and import to/from transmission system
    # ideally limited by transformer size (200 - 300 MVA would be reasonable for transmission transformer)
    # and perhaps trading regulations?
    @variables model begin
            import_export[i ∈ TRANSMISSION_NODES, t ∈ PERIODS]    ≥ 0
        0 ≤ export_to[i ∈ TRANSMISSION_NODES, t ∈ PERIODS]        ≤ 300
        0 ≤ import_from[i ∈ TRANSMISSION_NODES, t ∈ PERIODS]      ≤ 300
    end

    # import-export constraints
    # note: import_export is used to simplify nodal balance and objective function
    # when import_export ≥ 0, then the nodes are importing and imposed to trade cost (elprice)
    # otherwise, then the nodes are exporting and gain from trade cost (elprice)
    @constraint(model, 
                Trade[i ∈ TRANSMISSION_NODES, t ∈ PERIODS],
                import_export[i,t] == import_from[i,t] - export_to[i,t]
    )

    #=------------------------------------------------------------------------------
    VOLTAGE VARIABLES
    EQ (38)
    ------------------------------------------------------------------------------=#
    # voltage nominal and angles (in kV for voltage magnitude, degree for angle)
    # first nodes is excluded since it is considered slack bus
    @variables model begin
        0.9 * Vnom ≤ nodal_voltage[i ∈ NODES, t ∈ PERIODS] ≤ 1.1 * Vnom 
        nodal_angle[i ∈ NODES, t ∈ PERIODS]
    end

    # Slack bus voltage and angle over time, assumes the first entry as slack

    for t ∈ PERIODS
        fix(nodal_voltage[NODES[1], t], 
            Vnom,
            force=true
        )

        fix(nodal_angle[NODES[1], t], 
            0.0,
            force=true
        )

    end

    # @constraint(model, 
    #     Slack_Voltage[t ∈ PERIODS],
    #         nodal_voltage[NODES[1], t] == Vnom
    # )

    # @constraint(model, 
    #     Slack_Angle[t ∈ PERIODS],
    #         nodal_angle[NODES[1], t] == 0
    # )    

    #=------------------------------------------------------------------------------
    POWER FLOW VARIABLES
    EQ (40) - (41)
    ------------------------------------------------------------------------------=#
    # applied to lines and node from/to sets
    # since the flow can only be exist in power lines
    # but nodes still need to be taken into account
    # refer to Allard et el. (2020) eq. (4) - (7)
    # https://doi.org/10.1016/j.apenergy.2020.114958

    # apparent power limits
    @variables model begin
        -Lines_props[l][:s_max] ≤ active_flow[l ∈ LINES, t ∈ PERIODS]     ≤ Lines_props[l][:s_max]
        -Lines_props[l][:s_max] ≤ reactive_flow[l ∈ LINES, t ∈ PERIODS]   ≤ Lines_props[l][:s_max]
    end

    #=------------------------------------------------------------------------------
    Return
    ------------------------------------------------------------------------------=#

    Vars = (; 
            existing_generation,
            generation_investment, 
            storage_investment, 
            active_generation, 
            reactive_generation, 
            generation_spin,
            generation_on,
            gen_startup_cost,
            gen_partload_cost,
            gen_startup_CO2,
            gen_partload_CO2,            
            storage_charge, 
            storage_discharge, 
            storage_level, 
            nodal_voltage, 
            nodal_angle,
            import_export, 
            export_to, 
            import_from, 
            active_flow, 
            reactive_flow
    )

    return Vars

end         # end make_Variables


function make_Constraints(
    model::Model,
    sets::NamedTuple,
    params::NamedTuple,
    vars::NamedTuple,
    grid_infra::NamedTuple,
    profiles::NamedTuple
)

#=------------------------------------------------------------------------------
----------------------------- CONSTRAINTS --------------------------------------

Define Constraints of the model
current constraints:
1. Objective
2. power balance
3. storage balance
4. power flow
5. limits
6. ...

------------------------------------------------------------------------------=#
    ## Sets and parameters

    @unpack NODES, 
            TRANSMISSION_NODES,
            SE3_TRANS_NODES,
            NO1_TRANS_NODES,
            DK1_TRANS_NODES,
            COAST_NODES,
            GBG,
            GEN_TECHS, 
            EL_GEN, 
            HEAT_GEN, 
            H2_GEN, 
            STO_TECHS, 
            EL_STO, 
            HEAT_STO, 
            H2_STO, 
            PERIODS, 
            LINES, 
            NODE_FROM, 
            NODE_TO, 
            CHP,
            FC,
            WIND,
            PV,
            HP,
            BOILER,
            EC,
            FLEX_TH,
            THERMAL_1H,
            THERMAL_2H,
            # THERMAL_8H,
            THERMAL_12H = sets
    
    @unpack SE3_price,
            NO1_price,
            DK1_price, 
            Gentech_data, 
            Stotech_data, 
            Eldemand_data, 
            Heatdemand_data, 
            H2demand_data, 
            Discount_rate,
            Gen_cos_ϕ, 
            Gen_sin_ϕ, 
            Demand_cos_ϕ, 
            Demand_sin_ϕ, 
            Lines_props,
            Vnom = params

    ## Variables

    @unpack existing_generation,
            generation_investment, 
            storage_investment, 
            active_generation, 
            reactive_generation, 
            generation_spin,
            generation_on,
            gen_startup_cost,
            gen_partload_cost,
            gen_startup_CO2,
            gen_partload_CO2, 
            storage_charge, 
            storage_discharge, 
            storage_level, 
            nodal_voltage, 
            nodal_angle,
            import_export, 
            export_to, 
            import_from, 
            active_flow, 
            reactive_flow = vars

    #=------------------------------------------------------------------------------
    MODEL CONSTRAINTS
    ------------------------------------------------------------------------------=#

    #=------------------------------------------------------------------------------
    OBJECTIVES
    EQ (1)
    ------------------------------------------------------------------------------=#

    # Cost-related variables
    @variables model begin
        total_cost          # in k€
        capex               # in k€
        fix_om              # in k€
        fuel_cost           # in €
        var_om              # in €
        start_part_costs    # in €
        exp_imp_costs       # in €
        tax_cost            # in €
    end

    # define investments of technologies
    # capacity substracted by the lower_bound, which is the acquired existing data
    @expressions model begin
        CRF_gen[x ∈ GEN_TECHS],
            Discount_rate / (1 - 1/(1+Discount_rate)^Gentech_data[x].Lifetime)

        CRF_sto[s ∈ STO_TECHS],
            Discount_rate / (1 - 1/(1+Discount_rate)^Stotech_data[s].Lifetime)
    end

    # TODO: define more variables/expressions for costs
    # CAPEX costs (in k€/MW)
    @constraint(model, 
        capex == 
        # generation techs capacity investment costs
        sum(generation_investment[i, x] * Gentech_data[x].InvCost * CRF_gen[x] for i ∈ NODES, x ∈ GEN_TECHS) +
        # storage techs capacity investment costs
        sum(storage_investment[i, s] * Stotech_data[s].InvCost * CRF_sto[s] for i ∈ NODES, s ∈ STO_TECHS)                                             
    )

    # fix O&M costs (in k€/MW)
    @constraint(model, 
        fix_om ==
        # generation tech fix O&M costs, electrolyser defined differently since the fix OM is based on % of investment
        sum(generation_investment[i, x] * Gentech_data[x].FixOM for i ∈ NODES, x ∈ GEN_TECHS if x ∉ EC ) +                                              
        # generation tech fix O&M costs for electrolyser
        sum((generation_investment[i, x] * Gentech_data[x].InvCost * CRF_gen[x]) * Gentech_data[x].FixOM for i ∈ NODES, x ∈ GEN_TECHS if x ∈ EC ) +    
        # storage tech fix O&M costs, vanadium redox defined differently since the fix OM is based on % of investment
        sum(storage_investment[i, s] * Stotech_data[s].FixOM for i ∈ NODES, s ∈ STO_TECHS if s != "VRBAT") +                                            
        # storage tech fix O&M costs for vanadium redox
        sum((storage_investment[i, "VRBAT"] * Stotech_data["VRBAT"].InvCost * CRF_sto["VRBAT"]) * Stotech_data["VRBAT"].FixOM for i ∈ NODES)
    )

    # fuel costs (in €/MWh)
    @constraint(model, 
        fuel_cost == 
        # operational costs for techs that use fuel 
        sum( sum(active_generation[i, x, t] * Gentech_data[x].FuelPrice / Gentech_data[x].Efficiency for i ∈ NODES, x ∈ GEN_TECHS if x ∉ [HP, EC, "EB"] ) for t ∈ PERIODS) +
        # operational costs for HP and electrolyser uses elprice
        sum( sum(active_generation[i, x, t] * SE3_price[t] / Gentech_data[x].Efficiency for i ∈ NODES, x ∈ GEN_TECHS if x ∈ [HP, EC, "EB"] ) for t ∈ PERIODS)
    )

    # variable O/M costs (in €/MWh)
    @constraint(model, 
        var_om ==
        # generation tech variable O&M costs
        sum(active_generation[i, x, t] * Gentech_data[x].VarOM for i ∈ NODES, x ∈ GEN_TECHS, t ∈ PERIODS) +                                             
        # storage tech variable O&M costs 
        sum(storage_discharge[i, s, t] * Stotech_data[s].VarOM for i ∈ NODES, s ∈ STO_TECHS, t ∈ PERIODS)
    )

    # start-up/part load costs (in €/MWh)
    @constraint(model, 
        start_part_costs ==
        # startup costs
        sum(gen_startup_cost[i, x, t] for i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS) +
        # partload costs
        sum(gen_partload_cost[i, x, t] for i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS)
    )

    # import/export costs (in €/MWh)
    @constraint(model, 
        exp_imp_costs ==
        # export/import from transmission system
        sum(i ∈ SE3_TRANS_NODES ? import_export[i, t] * SE3_price[t] : 0 for i ∈ NODES, t ∈ PERIODS) +
        sum(i ∈ NO1_TRANS_NODES ? import_export[i, t] * NO1_price[t] : 0 for i ∈ NODES, t ∈ PERIODS) +
        sum(i ∈ DK1_TRANS_NODES ? import_export[i, t] * DK1_price[t] : 0 for i ∈ NODES, t ∈ PERIODS)
    )

    # taxes by using el for heat (in €/MWh)
    El_Heat_Tax = 60        # assumed tax in €/MWh for using el to generate heat
    @constraint(model, 
        tax_cost ==
        # taxes for EB
        sum(El_Heat_Tax * active_generation[i, "EB", t] / Gentech_data["EB"].Efficiency for i ∈ NODES, t ∈ PERIODS) +                                             
        # taxes for HP 
        sum(El_Heat_Tax * active_generation[i, x, t] / Gentech_data[x].Alpha for i ∈ NODES, x ∈ HP, t ∈ PERIODS)
    )

    # System Cost
    # represented in k€ for the total cost
    @constraint(model, SystemCost,
        total_cost * 1e3 == capex * 1e3 +       # in k€
                            fix_om * 1e3 +      # in k€
                            fuel_cost +         # in €
                            var_om +            # in €
                            start_part_costs +  # in €
                            exp_imp_costs +     # in €
                            tax_cost            # in €
    )

    #=------------------------------------------------------------------------------
    GENERATION LIMITS
    EQ (9) - (11)
    ------------------------------------------------------------------------------=#
    # general active and reactive generation limited by respective capacities
    # reactive generation is limited by the power factor relationship between active - reactive
    # reactive generation is assumed only applies for electricity generation technologies
    # Wind, PV, FLEX_TH are excluded because they are defined differently
    @constraints model begin
        # # active power generation
        # # commented because generation can opt not to generate
        # # perhaps any must run unit?
        # Active_Generation_Limit_lo[i ∈ NODES, x ∈ GEN_TECHS, t ∈ PERIODS; x != [WIND, PV, FLEX_TH]],
        #     existing_generation[i, x] ≤ active_generation[i, x, t]

        Active_Generation_Limit_up[i ∈ NODES, x ∈ GEN_TECHS, t ∈ PERIODS; x != [WIND, PV, FLEX_TH]],
            active_generation[i, x, t] ≤ existing_generation[i, x] + generation_investment[i, x]    
            
        # lower bound of reactive generation
        Reactive_Generation_Limit_lo[i ∈ NODES, x ∈ EL_GEN, t ∈ PERIODS],
            active_generation[i, x, t] * Gen_sin_ϕ[2] / Gen_cos_ϕ[2] ≤ reactive_generation[i, x, t]

        # upper bound of reactive generation
        Reactive_Generation_Limit_up[i ∈ NODES, x ∈ EL_GEN, t ∈ PERIODS],
            reactive_generation[i, x, t] ≤ active_generation[i, x, t] * Gen_sin_ϕ[1] / Gen_cos_ϕ[1]
    end

    # upper and lower bounds for generation capacity 
    # detailes are defined in set_gen_bounds function
    # include limits due to:
    # 1. existing power plants as lower bound
    # 2. voronoi area limits as upper limit
    # TODO: 3. RE potential
    set_gen_bounds(
        model,
        sets,
        params,
        vars,
        grid_infra
        )

    # specific generation limits
    # RE generations affected by the profile
    # based on Renewables Ninja profile
    @constraints model begin
        # offshore wind
        WOFF_Gen_Limitp[i ∈ COAST_NODES, x ∈ WIND, t ∈ PERIODS; x == "WOFF"],
            active_generation[i, x, t] ≤ profiles.WToff[t, i] * ( existing_generation[i, x] + generation_investment[i, x] )
            
        # onshore wind
        WON_Gen_Limit[i ∈ NODES, x ∈ WIND, t ∈ PERIODS; x == "WON"],
            active_generation[i, x, t] ≤ profiles.WTon[t, i] * ( existing_generation[i, x] + generation_investment[i, x] )

        # rooftop pv
        PV_roof_Limit[i ∈ NODES, x ∈ PV, t ∈ PERIODS; x == "PVROOF"],
            active_generation[i, x, t] ≤ profiles.PVfix[t, i] * ( existing_generation[i, x] + generation_investment[i, x] )

        # utility pv
        PV_util_Limit[i ∈ NODES, x ∈ PV, t ∈ PERIODS; x == "PVUTIL"],
            active_generation[i, x, t] ≤ profiles.PVfix[t, i] * ( existing_generation[i, x] + generation_investment[i, x] )

        # tracking pv
        PV_track_Limit[i ∈ NODES, x ∈ PV, t ∈ PERIODS; x == "PVTRACK"],
            active_generation[i, x, t] ≤ profiles.PVopt[t, i] * ( existing_generation[i, x] + generation_investment[i, x] )
    end

    # constraints related to two-variable approach
    # detailes are defined in flex_lim function
    flex_lim(
        model,
        sets,
        params,
        vars,
        )

    #=------------------------------------------------------------------------------
    STORAGE BALANCE LIMITS
    EQ (12) - (23)
    ------------------------------------------------------------------------------=#
    # storage-related constraints
    # Current initial storage level assumed to be 0
    initSto = zeros(length(NODES), length(STO_TECHS))#, length(PERIODS))
    Initial_Storage = AxisArrays.AxisArray(
                                    initSto, 
                                    AxisArrays.Axis{:node_id}(NODES), 
                                    AxisArrays.Axis{:Tech}(STO_TECHS)
    )#, Axis{:time}(PERIODS)) 

    @constraints model begin
        # Storage level limited by the capacity
        Storage_Level_Limit[i ∈ NODES, s ∈ STO_TECHS, t ∈ PERIODS],    
            storage_level[i, s, t] ≤ storage_investment[i, s]
    
        # Hourly storage level
        Storage_Balance[i ∈ NODES, s ∈ STO_TECHS, t ∈ PERIODS],
            storage_level[i, s, t] == 
            (t > 1 ? storage_level[i, s, t-1] : Initial_Storage[i,s]) + 
            storage_charge[i, s, t] * Stotech_data[s].Ch_eff - 
            storage_discharge[i, s, t] * Stotech_data[s].Dch_eff
    
        # Storage charge limited by the capacity and discharging rate
        Storage_charge_limit[i ∈ NODES, s ∈ STO_TECHS, t ∈ PERIODS],    
            storage_discharge[i, s, t] ≤ storage_investment[i, s] / Stotech_data[s].InjectionRate
    
        # Storage discharge limited by the capacity and charging rate
        Storage_discharge_limit[i ∈ NODES, s ∈ STO_TECHS, t ∈ PERIODS],    
            storage_discharge[i, s, t] ≤ storage_investment[i, s] / Stotech_data[s].WithdrawalRate    
        #TODO: cycle limits? for example line caverns?
        #TODO: losses in the storage? thermal, battery capacity, etc.?
    end

    #=------------------------------------------------------------------------------
    NODAL ENERGY BALANCES
    ------------------------------------------------------------------------------=#
    # Electricity nodal balance
    for node ∈ NODES, t ∈ PERIODS
        # Assuming the reactive demand corresponds to 0.95 cos phi
        Reactive_Demand = Eldemand_data .* Demand_sin_ϕ ./ Demand_cos_ϕ 

        # define electricity flow  from lines coming in/out of nodes
        # enter considered as generation/supply, vice versa
        # p denotes active power, q reactive power
        p_enter = sum(active_flow[line, t] for (idx, line) ∈ enumerate(LINES) if NODE_TO[idx] == node; init=0)
        p_exit = sum(active_flow[line, t] for (idx, line) ∈ enumerate(LINES) if NODE_FROM[idx] == node; init=0)
        q_enter = sum(reactive_flow[line, t] for (idx, line) ∈ enumerate(LINES) if NODE_TO[idx] == node; init=0)
        q_exit = sum(reactive_flow[line, t] for (idx, line) ∈ enumerate(LINES) if NODE_FROM[idx] == node; init=0)

        # EQ (31) #
        # active power nodal balance
        @constraint(model, 
            Eldemand_data[t, node] +                                                                                # el demand
            sum(active_generation[node, x, t] / Gentech_data[x].Alpha for x ∈ HP) +                            # for HP
            sum(active_generation[node, x, t] / Gentech_data[x].Efficiency for x ∈ BOILER) +                        # for boilers
            sum(active_generation[node, x, t] / Gentech_data[x].Efficiency for x ∈ EC) +                                                         # for electrolyser / H2 demand
            sum(storage_charge[node, s, t] for s ∈ EL_STO) +                                                        # charge battery
            sum(storage_charge[node, s, t] for s ∈ H2_STO) +                                                        # charge  h2 storage
            p_exit ≤                                                                                                # el flow to other nodes
            sum(active_generation[node, x, t] for x ∈ EL_GEN) +                        # el generation (active)
            sum(storage_discharge[node, s, t] for s ∈ EL_STO) +                                                     # battery discharge
            sum(active_generation[node, x, t] for x ∈ FC) +                                  # this applies for Fuel Cell (H2 -> EL)
            p_enter +                                                                                               # el flow to this node
            (node in TRANSMISSION_NODES ? import_export[node, t] : 0)                                               # import/export
        )
        
        # EQ (32) #
        # reactive power nodal balance     
        @constraint(model, 
            Reactive_Demand[t, node] +                                                  # reactive demand
            q_exit ≤                                                                    # reactive flow to other nodes
            sum(reactive_generation[node, x, t] for x ∈ EL_GEN) +                       # el generation (reactive)
            q_enter                                                                     # reactive flow to this node
        )
    end

    # EQ (28) - (29) #
    # Heat nodal balance
    for node ∈ NODES, t ∈ PERIODS
        # heat pipe flow equations would come here if there is any...
        @constraint(model,
            # efficiencies / el-heat conversion rate for HP and Boilers have been considered in el balance therefore not included
            Heatdemand_data[t, node] +
            sum(storage_charge[node, s, t] for s ∈ HEAT_STO) ≤                                 # charging heat storage
            sum(active_generation[node, x, t] / Gentech_data[x].Alpha for x ∈ CHP) +           # generation from heat techs, CHP if with alpha
            sum(active_generation[node, x, t] for x ∈ HP) +                                    # for HP
            sum(active_generation[node, x, t] for x ∈ BOILER) +                                # for boilers            
            sum(storage_discharge[node, s, t] for s ∈ HEAT_STO)                                # discharge from heat storage
            # possibility to buy heat from other region?
        )
    end

    # EQ (30) #
    # H2 nodal balance
    for node ∈ NODES, t ∈ PERIODS
        # H2 pipe flow equations would come here if there is any...
        @constraint(model,
            H2demand_data[t, node] +
            sum(storage_charge[node, s, t] for s ∈ H2_STO) +                                   # charging heat storage
            sum(active_generation[node, x, t] / Gentech_data[x].Efficiency for x ∈ FC) ≤       # fuel cell to convert h2 - el
            sum(active_generation[node, x, t] * Gentech_data[x].Efficiency for x ∈ EC) +       # electrolyser to convert el - h2
            sum(storage_discharge[node, s, t] for s ∈ H2_STO)                                  # discharge from H2 storage
            # possibility to buy H2 from other region?
        )
    end

    #=------------------------------------------------------------------------------
    LINEARISED POWER FLOW LIMITS
    ------------------------------------------------------------------------------=#
    # refer to Allard et el. (2020) eq. (8)
    # https://doi.org/10.1016/j.apenergy.2020.114958
    # apparent power limits have been defined in Variables - power flow variables
    # power flow equations
    for (idx, line) ∈ enumerate(LINES), t ∈ PERIODS
        G = Lines_props[line][:g_total]
        B = Lines_props[line][:b_total]

        # EQ (36) - (37) #
        # active flow
        @constraint(model,
            active_flow[line, t] == Vnom * (G * (nodal_voltage[NODE_FROM[idx], t] - nodal_voltage[NODE_TO[idx], t]) + 
                                    Vnom * B * (nodal_angle[NODE_TO[idx], t] - nodal_angle[NODE_FROM[idx], t]))
        )

        # reactive flow
        @constraint(model,
            reactive_flow[line, t] == Vnom * (B * (nodal_voltage[NODE_TO[idx], t] - nodal_voltage[NODE_FROM[idx], t]) + 
                                      Vnom * G * (nodal_angle[NODE_TO[idx], t] - nodal_angle[NODE_FROM[idx], t]))
        )

        # voltage angle difference limits (in degrees)
        # assumed to be limited by 30 deg
        @constraint(model,
            nodal_angle[NODE_FROM[idx], t] - nodal_angle[NODE_TO[idx], t] ≤ 30
        )

        @constraint(model,
            nodal_angle[NODE_FROM[idx], t] - nodal_angle[NODE_TO[idx], t] ≥ -30
        )

        # EQ (42) #
        # linearised thermal constraints
        @constraint(model,
            active_flow[line, t] + reactive_flow[line, t] ≤ sqrt(2) * Lines_props[line][:s_max]
        )

        @constraint(model,
            active_flow[line, t] - reactive_flow[line, t] ≤ sqrt(2) * Lines_props[line][:s_max]
        )

        @constraint(model,
            -active_flow[line, t] + reactive_flow[line, t] ≤ sqrt(2) * Lines_props[line][:s_max]
        )

        @constraint(model,
            -active_flow[line, t] - reactive_flow[line, t] ≤ sqrt(2) * Lines_props[line][:s_max]
        )
    end

    #=------------------------------------------------------------------------------
    CO2 LIMITS
    EQ (43)
    ------------------------------------------------------------------------------=#
    # CO2 limits and constraints in tonne
    # 2050 assumes net zero is achieved
    @constraint(model,
        sum(active_generation[i, x, t] * Gentech_data[x].Emission / Gentech_data[x].Efficiency for i ∈ NODES, x ∈ GEN_TECHS, t ∈ PERIODS) + 
        # sum(storage_discharge[i, s, t] * Stotech_data[s].Emission for i ∈ NODES, s ∈ STO_TECHS, t ∈ PERIODS) +
        sum(gen_startup_CO2[i, x, t] for i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS) +
        sum(gen_partload_CO2[i, x, t] for i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS) ≤
        0                                                                                      
    )

    Vars = (; 
            total_cost,
            existing_generation,
            generation_investment, 
            storage_investment, 
            active_generation, 
            reactive_generation, 
            generation_spin,
            generation_on,
            gen_startup_cost,
            gen_partload_cost,
            gen_startup_CO2,
            gen_partload_CO2,            
            storage_charge, 
            storage_discharge, 
            storage_level, 
            nodal_voltage, 
            nodal_angle,
            import_export, 
            export_to, 
            import_from, 
            active_flow, 
            reactive_flow
        )

    return Vars

end     # end make_Constraints


function set_gen_bounds(
    model::Model,
    sets::NamedTuple,
    params::NamedTuple,
    vars::NamedTuple,
    grid_infra::NamedTuple
)

    ## Sets and parameters

    @unpack NODES, 
            COAST_NODES,
            GBG, 
            GEN_TECHS, 
            STO_TECHS, 
            PERIODS = sets
    
    @unpack Gentech_data, 
            Stotech_data = params

    ## Variables

    @unpack existing_generation,
            generation_investment, 
            storage_investment, 
            active_generation = vars
    #=------------------------------------------------------------------------------
    GENERATION INVESTMENT LOWER BOUNDS

    based on the acquired power plant info from OpenStreetMap
    ------------------------------------------------------------------------------=#

    # define new columns that sums the el and heat existing capacity
    grid_infra.pp[!, "total_MW"] = grid_infra.pp.el_MW .+ grid_infra.pp.heat_MW

    # df only for the node, tech, and sum of el and heat MW
    aggregated_pp_df = combine(groupby(grid_infra.pp, [:node_id, :tech]), :total_MW => sum => :sum_MW)

    # empty dict and loop over the bounds    
    capacity_lower_bounds = Dict{Tuple{String, String}, Float64}()
    for row ∈ eachrow(aggregated_pp_df)
        capacity_lower_bounds[(row.node_id, row.tech)] = row.sum_MW
    end

    # define the lower bounds
    # set 0 as default value if no values are assigned
    for node ∈ NODES, tech ∈ GEN_TECHS
        fix(
            existing_generation[node, tech], 
            get(capacity_lower_bounds, (node, tech), 0.0),
            force=true
            )
    end

    # for node ∈ NODES, tech ∈ GEN_TECHS
    #     set_lower_bound(existing_generation[node, tech], get(capacity_lower_bounds, (node, tech), 0.0))
    # end

    #=------------------------------------------------------------------------------
    GENERATION INVESTMENT UPPER BOUNDS OR INELIGIBILITY
    ------------------------------------------------------------------------------=#

    # EQ (8) #
    # limits of area for each node based on Voronoi cell
    # assume 10% of the available area in 1000 m2
    # need to be better represented!!
    node_area = Dict(row.node_id => row.area_m2 for row ∈ eachrow(grid_infra.subs))

    # for node ∈ NODES
    #     @constraint(model,
    #         sum(( existing_generation[node, gen_tech] + generation_investment[node, gen_tech] ) * Gentech_data[gen_tech].Space_m2 for gen_tech ∈ GEN_TECHS) +
    #         sum(storage_investment[node, sto_tech] * Stotech_data[sto_tech].Space_m2 for sto_tech ∈ STO_TECHS) ≤ 
    #         node_area[node] #* 0.1
    #     )
    # end
    
    # nodes not in coastal municipalities not eligible for offshore wind farms
    for node ∈ NODES, t ∈ PERIODS
        if node ∉ COAST_NODES
            @constraint(model, 
                generation_investment[node, "WOFF"] == 0
            )

            @constraint(model, 
                active_generation[node, "WOFF", t] == 0
            )
        end
    end

    # techs that most probably not invested in future
    # for now still include waste CHP (WCHP) and waste boiler (WBO)
    TECHS_NOT_INVESTED = [
                        "COCHP", 
                        # "WCHP", 
                        # "WBO", 
                        "HYD"
    ]

    # further limitations for coal, waste CHP, waste boiler
    # as they are not likely be further investment in future
    # set existing value as upper bound if any
    # otherwise set as 0
    for node ∈ NODES, tech ∈ TECHS_NOT_INVESTED
        fix(generation_investment[node, tech], 
            0.0, 
            force=true
        )
    end

    
    # Pit Thermal Storage is not feasible in Gothenburg
    for node ∈ GBG
        fix(storage_investment[node, "PTES"], 
            0.0,
            force=true
        )
    end

    # Seawater Heat Pumps could only be invested in coastal area
    # does not make sense to buy sea water and transport it
    for node ∈ NODES
        if node ∉ COAST_NODES
            fix(generation_investment[node, "HPSW"], 
                0.0, 
                force=true
            )
        end
    end

end         # end set_gen_bounds


function flex_lim(
    model::Model,
    sets::NamedTuple,
    params::NamedTuple,
    vars::NamedTuple,
)
#=------------------------------------------------------------------------------
-------------------------- FLEX LIM CONSTRAINTS --------------------------------

Define Constraints of the model
related to the two-variable approach:
1. gen spin
2. start up cost
3. part load cost
4. ...
------------------------------------------------------------------------------=#
    ## Sets and parameters

    @unpack NODES, 
            GEN_TECHS, 
            PERIODS, 
            FLEX_TH,
            THERMAL_1H,
            THERMAL_2H,
            # THERMAL_8H,
            THERMAL_12H = sets
    
    @unpack Gentech_data = params

    ## Variables

    @unpack existing_generation,
            generation_investment, 
            active_generation, 
            generation_spin,
            generation_on,
            gen_startup_cost,
            gen_partload_cost,
            gen_startup_CO2,
            gen_partload_CO2 = vars

    #=------------------------------------------------------------------------------
    TWO-VARIABLE APPROACH CONSTRAINTS
    EQ (24) - (27)
    ------------------------------------------------------------------------------=#
    # maximum gen bounded by available "hot capacity" (spin)
    @constraint(model, Spin_lim[i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS],
        generation_spin[i, x, t] ≤ ( existing_generation[i, x] + generation_investment[i, x] )
    )
    
    @constraint(model, Ramping_up[i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS],
        active_generation[i, x, t] ≤ generation_spin[i, x, t]
    )

    # minimum gen bounded by minimum load
    @constraint(model, Ramping_down[i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS],
        Gentech_data[x].MinLoad * generation_spin[i, x, t] ≤ active_generation[i, x, t]
    )

    # start-up limits
    @constraint(model, [i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS; t ≥ 2],
        generation_on[i, x, t] ≥ 
        generation_spin[i, x, t] - generation_spin[i, x, t-1]
    )

    # startup cost and the CO2 emission at start up
    @constraint(model, Startup_cost[i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS],
        gen_startup_cost[i, x, t] ≥ Gentech_data[x].StartCost * generation_on[i, x, t]
    )

    @constraint(model, Startup_CO2[i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS],
        gen_startup_CO2[i, x, t] ≥ Gentech_data[x].StartCO2 * generation_on[i, x, t]
    )

    # part load cost and the CO2 emission by part load
    @constraint(model, Partload_cost[i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS],
        gen_partload_cost[i, x, t] ≥ 
        Gentech_data[x].PartLoadCost * (generation_spin[i, x, t] - active_generation[i, x, t])
    )

    @constraint(model, Partload_CO2[i ∈ NODES, x ∈ FLEX_TH, t ∈ PERIODS],
        gen_partload_CO2[i, x, t] ≥ 
        Gentech_data[x].PartLoadCO2 * (generation_spin[i, x, t] - active_generation[i, x, t])
    )

    # Ramping limits
    # should apply when the period is longer than a day
    if size(PERIODS, 1) ≥ 12
        # new tech sets based on the start up time duration
        # defined manually based on data in excel

        @constraint(model, Ramp_1h[i ∈ NODES, x ∈ THERMAL_1H, t ∈ PERIODS; t ≥ 2],
            generation_on[i, x, t] ≤ 
            ( existing_generation[i, x] + generation_investment[i, x] ) - generation_spin[i, x, t-1]
        )

        @constraint(model, Ramp_2h[i ∈ NODES, x ∈ THERMAL_2H, t ∈ PERIODS; t ≥ 3],
            generation_on[i, x, t] ≤ 
            ( existing_generation[i, x] + generation_investment[i, x] ) - generation_spin[i, x, t-2]
        )

        # @constraint(model, Ramp_8h[i ∈ NODES, x ∈ THERMAL_8H, t ∈ PERIODS; t ≥ 9],
        #     generation_on[i, x, t] ≤ 
        #     ( existing_generation[i, x] + generation_investment[i, x] ) - generation_spin[i, x, t-8]
        # )

        @constraint(model, Ramp_12h[i ∈ NODES, x ∈ THERMAL_12H, t ∈ PERIODS; t ≥ 13],
            generation_on[i, x, t] ≤ 
            ( existing_generation[i, x] + generation_investment[i, x] ) - generation_spin[i, x, t-12]
        )
    end

end         # end flex_lim