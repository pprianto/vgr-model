function cost_constraints(
    model::Model,
    sets::ModelSets,
    params::ModelParameters,
    vars::ModelVariables,
    grid_infra::GridInfrastructures,
    profiles::Profiles
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
            Reactive_Demand, 
            Heatdemand_data, 
            H2demand_data, 
            # Discount_rate,
            Gen_cos_ϕ, 
            Gen_sin_ϕ, 
            Demand_cos_ϕ, 
            Demand_sin_ϕ, 
            Lines_props = params
            # Vnom = params

    ## Variables

    @unpack total_cost,
            capex,
            fix_om,
            fuel_cost,
            var_om,
            start_part_costs,
            exp_imp_costs,
            tax_cost,
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
            reactive_flow = vars

    #=------------------------------------------------------------------------------
    MODEL CONSTRAINTS
    ------------------------------------------------------------------------------=#

    #=------------------------------------------------------------------------------
    OBJECTIVES
    EQ (1)
    ------------------------------------------------------------------------------=#

    println("---------------------------------------------")
    println("Costs constraints")
    println("---------------------------------------------")

    @objective(model, Min, total_cost)

    # define investments of technologies
    # capacity substracted by the lower_bound, which is the acquired existing data
    @expressions model begin
        CRF_gen[x ∈ GEN_TECHS],
            options.Discount_rate / (1 - 1/(1+options.Discount_rate)^Gentech_data[x].Lifetime)

        CRF_sto[s ∈ STO_TECHS],
            options.Discount_rate / (1 - 1/(1+options.Discount_rate)^Stotech_data[s].Lifetime)
    end

    # CAPEX costs (in k€/MW)
    @constraint(model, 
        capex ≥ 
        # generation techs capacity investment costs
        sum( 
            sum(generation_investment[i, x] * Gentech_data[x].InvCost * CRF_gen[x] for x ∈ GEN_TECHS) +
        # storage techs capacity investment costs
            sum(storage_investment[i, s] * Stotech_data[s].InvCost * CRF_sto[s] for s ∈ STO_TECHS) 
        for i ∈ NODES  )                                      
    )

    # fix O&M costs (in k€/MW)
    @constraint(model, 
        fix_om ≥
        # generation tech fix O&M costs, electrolyser defined differently since the fix OM is based on % of investment
        sum(
            sum(generation_investment[i, x] * Gentech_data[x].FixOM for x ∈ GEN_TECHS if x ∉ EC ) +                                              
            # generation tech fix O&M costs for electrolyser
            sum((generation_investment[i, x] * Gentech_data[x].InvCost * CRF_gen[x]) * Gentech_data[x].FixOM for x ∈ GEN_TECHS if x ∈ EC ) +    
            # storage tech fix O&M costs, vanadium redox defined differently since the fix OM is based on % of investment
            sum(storage_investment[i, s] * Stotech_data[s].FixOM for s ∈ STO_TECHS if s != "VRBAT") +                                            
            # storage tech fix O&M costs for vanadium redox
            sum((storage_investment[i, "VRBAT"] * Stotech_data["VRBAT"].InvCost * CRF_sto["VRBAT"]) * Stotech_data["VRBAT"].FixOM)
        for i ∈ NODES ) 
    )

    # fuel costs (in €/MWh)
    @constraint(model, 
        fuel_cost ≥ 
        # operational costs for techs that use fuel 
        sum( 
            sum(
                sum(active_generation[t, i, x] * Gentech_data[x].FuelPrice / Gentech_data[x].Efficiency for x ∈ GEN_TECHS if x ∉ [HP, EC, "EB"]) +
                # operational costs for HP and electrolyser uses elprice
                sum(active_generation[t, i, x] * SE3_price[t] / Gentech_data[x].Efficiency for x ∈ GEN_TECHS if x ∈ [HP, EC, "EB"] )
                for i ∈ NODES)
        for t ∈ PERIODS)
    )

    # variable O/M costs (in €/MWh)
    @constraint(model, 
        var_om ≥
        # generation tech variable O&M costs
        sum( 
            sum( 
                sum(active_generation[t, i, x] * Gentech_data[x].VarOM for x ∈ GEN_TECHS) +                                             
                # storage tech variable O&M costs 
                sum(storage_discharge[t, i, s] * Stotech_data[s].VarOM for s ∈ STO_TECHS)
                for i ∈ NODES)
        for t ∈ PERIODS)
    )

    # start-up/part load costs (in €/MWh)
    if options.FlexLim == :yes
        @constraint(model, 
            start_part_costs ≥
            # startup costs
            sum(
                sum(
                    sum(gen_startup_cost[t, i, x] for x ∈ FLEX_TH) +
                    # partload costs
                    sum(gen_partload_cost[t, i, x] for x ∈ FLEX_TH)
                for i ∈ NODES)
            for t ∈ PERIODS)
        )
    end

    # import/export costs (in €/MWh)
    @constraint(model, 
        exp_imp_costs ≥
        # export/import from transmission system
        sum( 
            sum(i ∈ SE3_TRANS_NODES ? import_export[t, i] * SE3_price[t] : 0 for i ∈ NODES) +
            sum(i ∈ NO1_TRANS_NODES ? import_export[t, i] * NO1_price[t] : 0 for i ∈ NODES) +
            sum(i ∈ DK1_TRANS_NODES ? import_export[t, i] * DK1_price[t] : 0 for i ∈ NODES) 
        for t ∈ PERIODS)
    )

    # taxes by using el for heat (in €/MWh)
    @constraint(model, 
        tax_cost ≥
        # taxes for EB and HP
        sum( 
            sum( 
                options.El_Heat_Tax * active_generation[t, i, "EB"] / Gentech_data["EB"].Efficiency +                                             
                # taxes for HP 
                sum(options.El_Heat_Tax * active_generation[t, i, x] / Gentech_data[x].Alpha for x ∈ HP) 
            for i ∈ NODES)
        for t ∈ PERIODS)
    )

    # System Cost
    # represented in k€ for the total cost
    if options.FlexLim == :yes
        @constraint(model, SystemCost,
            total_cost * 1e3 ≥ capex * 1e3 +       # in k€
                                fix_om * 1e3 +      # in k€
                                fuel_cost +         # in €
                                var_om +            # in €
                                start_part_costs +  # in €
                                exp_imp_costs +     # in €
                                tax_cost            # in €
        )

    else
        @constraint(model, SystemCost,
            total_cost * 1e3 ≥ capex * 1e3 +       # in k€
                                fix_om * 1e3 +      # in k€
                                fuel_cost +         # in €
                                var_om +            # in €
                                exp_imp_costs +     # in €
                                tax_cost            # in €
        )
    end

    #=------------------------------------------------------------------------------
    CO2 LIMITS
    EQ (43)
    ------------------------------------------------------------------------------=#
    # CO2 limits and constraints in tonne
    # 2050 assumes net zero is achieved
    println("---------------------------------------------")
    println("CO2 constraints")
    println("---------------------------------------------")

    if options.FlexLim == :yes
        @constraint(model,
            sum( 
                sum( 
                    sum(active_generation[t, i, x] * Gentech_data[x].Emission / Gentech_data[x].Efficiency for x ∈ GEN_TECHS) +
                    sum(gen_startup_CO2[t, i, x] + gen_partload_CO2[t, i, x] for x ∈ FLEX_TH)
                for i ∈ NODES)
            for t ∈ PERIODS) ≤
            options.CO2_limit                                                                                      
        )
    else
        @constraint(model,
            sum( 
                sum( 
                    sum(active_generation[t, i, x] * Gentech_data[x].Emission / Gentech_data[x].Efficiency for x ∈ GEN_TECHS)
                for i ∈ NODES)
            for t ∈ PERIODS) ≤
            options.CO2_limit                                                                                      
        )
    end

end     # end cost_constraints


function gen_constraints(
    model::Model,
    sets::ModelSets,
    params::ModelParameters,
    vars::ModelVariables,
    grid_infra::GridInfrastructures,
    profiles::Profiles
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
            COAST_NODES,
            GBG,
            GEN_TECHS,
            EL_GEN, 
            PERIODS, 
            CHP,
            WIND,
            PV,
            FLEX_TH = sets
    
    @unpack Gentech_data, 
            Stotech_data, 
            Gen_cos_ϕ, 
            Gen_sin_ϕ = params

    ## Variables

    @unpack existing_generation,
            generation_investment,
            storage_investment,
            active_generation,
            reactive_generation = vars

    #=------------------------------------------------------------------------------
    MODEL CONSTRAINTS
    ------------------------------------------------------------------------------=#
    # upper and lower bounds for generation capacity 
    # detailes are defined in set_gen_bounds function
    # include limits due to:
    # 1. existing power plants as lower bound
    # 2. voronoi area limits as upper limit
    # TODO: 3. RE potential

    println("---------------------------------------------")
    println("Gen bounds constraints")
    println("---------------------------------------------")

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
    # for node ∈ NODES 
    #     for tech ∈ GEN_TECHS
    #         set_lower_bound(
    #             existing_generation[node, tech], 
    #             get(capacity_lower_bounds, (node, tech), 0.0)
    #             )
    #     end
    # end

    # since it is existing, set as fixed
    for node ∈ NODES 
        for tech ∈ GEN_TECHS
            fix(
                existing_generation[node, tech], 
                get(capacity_lower_bounds, (node, tech), 0.0),
                force=true
                )
        end
    end

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
    
    # Nodes not in coastal municipalities not eligible for offshore wind farms (WOFF)
    # Seawater Heat Pumps (HPSW) could only be invested in coastal area
    # does not make sense to buy sea water and transport it
    
    for t ∈ PERIODS
        for node ∈ COAST_NODES
            fix(
            active_generation[t, node, "WOFF"], 
            0.0, 
            force=true
            )
        end
    end

    for node ∈ COAST_NODES
        fix(
        generation_investment[node, "WOFF"], 
        0.0, 
        force=true
        )
        
        fix(generation_investment[node, "HPSW"], 
        0.0, 
        force=true
        )
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
    for node ∈ NODES
        for tech ∈ TECHS_NOT_INVESTED
            fix(generation_investment[node, tech], 
                0.0, 
                force=true
            )
        end
    end
    
    # Pit Thermal Storage is not feasible in Gothenburg
    for node ∈ GBG
        fix(storage_investment[node, "PTES"], 
            0.0,
            force=true
        )
    end


    # specific generation limits
    # RE generations affected by the profile
    # based on Renewables Ninja profile
    @constraints model begin
        # offshore wind
        WOFF_Gen_Limit[t ∈ PERIODS, i ∈ COAST_NODES],
            active_generation[t, i, "WOFF"] ≤ profiles.WToff[t, i] * ( existing_generation[i, "WOFF"] + generation_investment[i, "WOFF"] )
            
        # onshore wind
        WON_Gen_Limit[t ∈ PERIODS, i ∈ NODES],
            active_generation[t, i, "WON"] ≤ profiles.WTon[t, i] * ( existing_generation[i, "WON"] + generation_investment[i, "WON"] )

        # rooftop pv
        PV_roof_Limit[t ∈ PERIODS, i ∈ NODES],
            active_generation[t, i, "PVROOF"] ≤ profiles.PVfix[t, i] * ( existing_generation[i, "PVROOF"] + generation_investment[i, "PVROOF"] )

        # utility pv
        PV_util_Limit[t ∈ PERIODS, i ∈ NODES],
            active_generation[t, i, "PVUTIL"] ≤ profiles.PVfix[t, i] * ( existing_generation[i, "PVUTIL"] + generation_investment[i, "PVUTIL"] )

        # tracking pv
        PV_track_Limit[t ∈ PERIODS, i ∈ NODES],
            active_generation[t, i, "PVTRACK"] ≤ profiles.PVopt[t, i] * ( existing_generation[i, "PVTRACK"] + generation_investment[i, "PVTRACK"] )
    end

    println("---------------------------------------------")
    println("Gen limits constraints")
    println("---------------------------------------------")

    #=------------------------------------------------------------------------------
    GENERATION LIMITS
    EQ (9) - (11)
    ------------------------------------------------------------------------------=#
    # general active and reactive generation limited by respective capacities
    # reactive generation is limited by the power factor relationship between active - reactive
    # reactive generation is assumed only applies for electricity generation technologies
    # Wind, PV, FLEX_TH are excluded because they are defined differently
    
    if options.FlexLim == :yes
        @constraint(model, Active_Generation_Limit_up[t ∈ PERIODS, i ∈ NODES, x ∈ GEN_TECHS; x != [WIND, PV, FLEX_TH]],
            active_generation[t, i, x] ≤ existing_generation[i, x] + generation_investment[i, x])
    else
        @constraint(model, Active_Generation_Limit_up[t ∈ PERIODS, i ∈ NODES, x ∈ GEN_TECHS; x != [WIND, PV]],
        active_generation[t, i, x] ≤ existing_generation[i, x] + generation_investment[i, x])
    end  
    
    @constraints model begin
        # # active power generation
        # # commented because generation can opt not to generate
        # # perhaps any must run unit?
        # Active_Generation_Limit_lo[i ∈ NODES, x ∈ GEN_TECHS, t ∈ PERIODS; x != [WIND, PV, FLEX_TH]],
        #     existing_generation[i, x] ≤ active_generation[t, i, x]
        
        # lower bound of reactive generation
        Reactive_Generation_Limit_lo[t ∈ PERIODS, i ∈ NODES, x ∈ EL_GEN],
            active_generation[t, i, x] * Gen_sin_ϕ[2] / Gen_cos_ϕ[2] ≤ reactive_generation[t, i, x]

        # upper bound of reactive generation
        Reactive_Generation_Limit_up[t ∈ PERIODS, i ∈ NODES, x ∈ EL_GEN],
            reactive_generation[t, i, x] ≤ active_generation[t, i, x] * Gen_sin_ϕ[1] / Gen_cos_ϕ[1]
    end


end     # end gen_Constraints


function flex_lim(
    model::Model,
    sets::ModelSets,
    params::ModelParameters,
    vars::ModelVariables,
    grid_infra::GridInfrastructures,
    profiles::Profiles
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
            storage_investment,
            active_generation,
            generation_spin,
            generation_on,
            gen_startup_cost,
            gen_partload_cost,
            gen_startup_CO2,
            gen_partload_CO2 = vars

    #=------------------------------------------------------------------------------
    MODEL CONSTRAINTS
    ------------------------------------------------------------------------------=#
    # constraints related to two-variable approach
    # detailes are defined in flex_lim function

    println("---------------------------------------------")
    println("Flex lim constraints")
    println("---------------------------------------------")

    #=------------------------------------------------------------------------------
    TWO-VARIABLE APPROACH CONSTRAINTS
    EQ (24) - (27)
    ------------------------------------------------------------------------------=#
    # maximum gen bounded by available "hot capacity" (spin)
    @constraint(model, Spin_lim[t ∈ PERIODS, i ∈ NODES, x ∈ FLEX_TH],
        generation_spin[t, i, x] ≤ ( existing_generation[i, x] + generation_investment[i, x] )
    )
    
    @constraint(model, Ramping_up[t ∈ PERIODS, i ∈ NODES, x ∈ FLEX_TH],
        active_generation[t, i, x] ≤ generation_spin[t, i, x]
    )

    # minimum gen bounded by minimum load
    @constraint(model, Ramping_down[t ∈ PERIODS, i ∈ NODES, x ∈ FLEX_TH],
        Gentech_data[x].MinLoad * generation_spin[t, i, x] ≤ active_generation[t, i, x]
    )

    # start-up limits
    @constraint(model, Hourly_increment[t ∈ PERIODS, i ∈ NODES, x ∈ FLEX_TH; t ≥ 2],
        generation_on[t, i, x] ≥ 
        generation_spin[t, i, x] - generation_spin[t-1, i, x]
    )

    # startup cost and the CO2 emission at start up
    @constraint(model, Startup_cost[t ∈ PERIODS, i ∈ NODES, x ∈ FLEX_TH],
        gen_startup_cost[t, i, x] ≥ Gentech_data[x].StartCost * generation_on[t, i, x]
    )

    @constraint(model, Startup_CO2[t ∈ PERIODS, i ∈ NODES, x ∈ FLEX_TH],
        gen_startup_CO2[t, i, x] ≥ Gentech_data[x].StartCO2 * generation_on[t, i, x]
    )

    # part load cost and the CO2 emission by part load
    @constraint(model, Partload_cost[t ∈ PERIODS, i ∈ NODES, x ∈ FLEX_TH],
        gen_partload_cost[t, i, x] ≥ 
        Gentech_data[x].PartLoadCost * (generation_spin[t, i, x] - active_generation[t, i, x])
    )

    @constraint(model, Partload_CO2[t ∈ PERIODS, i ∈ NODES, x ∈ FLEX_TH],
        gen_partload_CO2[t, i, x] ≥ 
        Gentech_data[x].PartLoadCO2 * (generation_spin[t, i, x] - active_generation[t, i, x])
    )

    # Ramping limits
    # should apply when the period is longer than a day
    if size(PERIODS, 1) ≥ 12
        # new tech sets based on the start up time duration
        # defined manually based on data in excel

        @constraint(model, Ramp_1h[t ∈ PERIODS, i ∈ NODES, x ∈ THERMAL_1H; t ≥ 2],
            generation_on[t, i, x] ≤ 
            ( existing_generation[i, x] + generation_investment[i, x] ) - generation_spin[t-1, i, x]
        )

        @constraint(model, Ramp_2h[t ∈ PERIODS, i ∈ NODES, x ∈ THERMAL_2H; t ≥ 3],
            generation_on[t, i, x] ≤ 
            ( existing_generation[i, x] + generation_investment[i, x] ) - generation_spin[t-2, i, x]
        )

        # @constraint(model, Ramp_8h[t ∈ PERIODS, i ∈ NODES, x ∈ THERMAL_8H; t ≥ 9],
        #     generation_on[t, i, x] ≤ 
        #     ( existing_generation[i, x] + generation_investment[i, x] ) - generation_spin[i, x, t-8]
        # )

        @constraint(model, Ramp_12h[t ∈ PERIODS, i ∈ NODES, x ∈ THERMAL_12H; t ≥ 13],
            generation_on[t, i, x] ≤ 
            ( existing_generation[i, x] + generation_investment[i, x] ) - generation_spin[t-12, i, x]
        )
    end

end     # end flex_lim


function sto_constraints(
    model::Model,
    sets::ModelSets,
    params::ModelParameters,
    vars::ModelVariables,
    grid_infra::GridInfrastructures,
    profiles::Profiles
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
            STO_TECHS, 
            EL_STO, 
            HEAT_STO, 
            H2_STO, 
            PERIODS = sets
    
    @unpack Stotech_data = params
            # Vnom = params

    ## Variables

    @unpack storage_investment,
            storage_charge,
            storage_discharge,
            storage_level = vars

    #=------------------------------------------------------------------------------
    MODEL CONSTRAINTS
    ------------------------------------------------------------------------------=#
    #=------------------------------------------------------------------------------
    STORAGE BALANCE LIMITS
    EQ (12) - (23)
    ------------------------------------------------------------------------------=#
    # storage-related constraints
    # Current initial storage level assumed to be 0
    println("---------------------------------------------")
    println("Storage constraints")
    println("---------------------------------------------")

    initSto = zeros(length(NODES), length(STO_TECHS))#, length(PERIODS))
    Initial_Storage = AxisArrays.AxisArray(
                                    initSto, 
                                    AxisArrays.Axis{:node_id}(NODES), 
                                    AxisArrays.Axis{:Tech}(STO_TECHS)
    )#, Axis{:time}(PERIODS)) 

    @constraints model begin
        # Storage level limited by the capacity
        Storage_Level_Limit[t ∈ PERIODS, i ∈ NODES, s ∈ STO_TECHS],    
            storage_level[t, i, s] ≤ storage_investment[i, s]
    
        # Hourly storage level
        Storage_Balance[t ∈ PERIODS, i ∈ NODES, s ∈ STO_TECHS],
            storage_level[t, i, s] == 
            (t > 1 ? storage_level[t-1, i, s] : Initial_Storage[i,s]) + 
            storage_charge[t, i, s] * Stotech_data[s].Ch_eff - 
            storage_discharge[t, i, s] / Stotech_data[s].Dch_eff
    
        # Storage charge limited by the capacity and discharging rate
        Storage_charge_limit[t ∈ PERIODS, i ∈ NODES, s ∈ STO_TECHS],    
            storage_charge[t, i, s] ≤ storage_investment[i, s] / Stotech_data[s].InjectionRate
    
        # Storage discharge limited by the capacity and charging rate
        Storage_discharge_limit[t ∈ PERIODS, i ∈ NODES, s ∈ STO_TECHS],    
            storage_discharge[t, i, s] ≤ storage_investment[i, s] / Stotech_data[s].WithdrawalRate    
        
        # Line Rock cavern cycle limits
        Line_caverns_limit_1[i ∈ NODES],
            sum( storage_charge[t, i, "LRC"] * Stotech_data["LRC"].Ch_eff for t ∈ PERIODS) ≤ storage_investment[i, "LRC"] * 20

        Line_caverns_limit_2[i ∈ NODES],
            sum( storage_discharge[t, i, "LRC"] / Stotech_data["LRC"].Dch_eff for t ∈ PERIODS) ≤ storage_investment[i, "LRC"] * 20    

        #TODO: losses in the storage? thermal, battery capacity, etc.?
    end

end     # end sto_Constraints


function enbal_constraints(
    model::Model,
    sets::ModelSets,
    params::ModelParameters,
    vars::ModelVariables,
    grid_infra::GridInfrastructures,
    profiles::Profiles
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
            EC = sets
    
    @unpack Gentech_data, 
            Stotech_data, 
            Eldemand_data,
            Reactive_Demand, 
            Heatdemand_data, 
            H2demand_data = params

    ## Variables

    @unpack active_generation,
            reactive_generation,
            storage_charge,
            storage_discharge,
            storage_level,
            import_export,
            export_to,
            import_from,
            active_flow,
            reactive_flow = vars

    #=------------------------------------------------------------------------------
    MODEL CONSTRAINTS
    ------------------------------------------------------------------------------=#

    #=------------------------------------------------------------------------------
    NODAL ENERGY BALANCES
    ------------------------------------------------------------------------------=#
    println("---------------------------------------------")
    println("El nodal balance constraints")
    println("---------------------------------------------")

    # Electricity nodal balance
    for t ∈ PERIODS
        for node ∈ NODES

            # define electricity flow  from lines coming in/out of nodes
            # enter considered as generation/supply, vice versa
            # p denotes active power, q reactive power
            p_enter = sum(active_flow[t, line] for (idx, line) ∈ enumerate(LINES) if NODE_TO[idx] == node; init=0)
            p_exit = sum(active_flow[t, line] for (idx, line) ∈ enumerate(LINES) if NODE_FROM[idx] == node; init=0)
            q_enter = sum(reactive_flow[t, line] for (idx, line) ∈ enumerate(LINES) if NODE_TO[idx] == node; init=0)
            q_exit = sum(reactive_flow[t, line] for (idx, line) ∈ enumerate(LINES) if NODE_FROM[idx] == node; init=0)

            # EQ (31) #
            # active power nodal balance
            @constraint(model, 
                Eldemand_data[t, node] +                                                                                # el demand
                sum(active_generation[t, node, x] / Gentech_data[x].Alpha for x ∈ HP) +                            # for HP
                sum(active_generation[t, node, x] / Gentech_data[x].Efficiency for x ∈ BOILER) +                        # for boilers
                sum(active_generation[t, node, x] / Gentech_data[x].Efficiency for x ∈ EC) +                                                         # for electrolyser / H2 demand
                sum(storage_charge[t, node, s] for s ∈ EL_STO) +                                                        # charge battery
                sum(storage_charge[t, node, s] for s ∈ H2_STO) +                                                        # charge  h2 storage
                p_exit ≤                                                                                                # el flow to other nodes
                sum(active_generation[t, node, x] for x ∈ EL_GEN) +                        # el generation (active)
                sum(storage_discharge[t, node, s] for s ∈ EL_STO) +                                                     # battery discharge
                sum(active_generation[t, node, x] for x ∈ FC) +                                  # this applies for Fuel Cell (H2 -> EL)
                p_enter +                                                                                               # el flow to this node
                (node in TRANSMISSION_NODES ? import_export[t, node] : 0)                                               # import/export
            )
            
            # EQ (32) #
            # reactive power nodal balance     
            @constraint(model, 
                Reactive_Demand[t, node] +                                                  # reactive demand
                q_exit ≤                                                                    # reactive flow to other nodes
                sum(reactive_generation[t, node, x] for x ∈ EL_GEN) +                       # el generation (reactive)
                q_enter                                                                     # reactive flow to this node
            )
        end
    end

    println("---------------------------------------------")
    println("Heat nodal balance constraints")
    println("---------------------------------------------")

    # EQ (28) - (29) #
    # Heat nodal balance
    for t ∈ PERIODS
        for node ∈ NODES
            # heat pipe flow equations would come here if there is any...
            @constraint(model,
                # efficiencies / el-heat conversion rate for HP and Boilers have been considered in el balance therefore not included
                Heatdemand_data[t, node] +
                sum(storage_charge[t, node, s] for s ∈ HEAT_STO) ≤                                 # charging heat storage
                sum(active_generation[t, node, x] / Gentech_data[x].Alpha for x ∈ CHP) +           # generation from heat techs, CHP if with alpha
                sum(active_generation[t, node, x] for x ∈ HP) +                                    # for HP
                sum(active_generation[t, node, x] for x ∈ BOILER) +                                # for boilers            
                sum(storage_discharge[t, node, s] for s ∈ HEAT_STO)                                # discharge from heat storage
                # possibility to buy heat from other region?
            )
        end
    end

    println("---------------------------------------------")
    println("H2 nodal balance constraints")
    println("---------------------------------------------")

    # EQ (30) #
    # H2 nodal balance
    for t ∈ PERIODS
        for node ∈ NODES
            # H2 pipe flow equations would come here if there is any...
            @constraint(model,
                H2demand_data[t, node] +
                sum(storage_charge[t, node, s] for s ∈ H2_STO) +                                   # charging heat storage
                sum(active_generation[t, node, x] / Gentech_data[x].Efficiency for x ∈ FC) ≤       # fuel cell to convert h2 - el
                sum(active_generation[t, node, x] * Gentech_data[x].Efficiency for x ∈ EC) +       # electrolyser to convert el - h2
                sum(storage_discharge[t, node, s] for s ∈ H2_STO)                                  # discharge from H2 storage
                # possibility to buy H2 from other region?
            )
        end
    end

end     # end enbal_Constraints



function power_flow_constraints(
    model::Model,
    sets::ModelSets,
    params::ModelParameters,
    vars::ModelVariables,
    grid_infra::GridInfrastructures,
    profiles::Profiles
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
            PERIODS, 
            LINES, 
            NODE_FROM, 
            NODE_TO = sets
    
    @unpack Lines_props = params
            # Vnom = params

    ## Variables

    @unpack nodal_voltage,
            nodal_angle,
            active_flow,
            reactive_flow = vars

    #=------------------------------------------------------------------------------
    MODEL CONSTRAINTS
    ------------------------------------------------------------------------------=#
    #=------------------------------------------------------------------------------
    LINEARISED POWER FLOW LIMITS
    ------------------------------------------------------------------------------=#
    # refer to Allard et el. (2020) eq. (8)
    # https://doi.org/10.1016/j.apenergy.2020.114958
    # apparent power limits have been defined in Variables - power flow variables
    # power flow equations

    println("---------------------------------------------")
    println("Power flow constraints")
    println("---------------------------------------------")

    if options.power_flow == :lac

        for t ∈ PERIODS 
            for (idx, line) ∈ enumerate(LINES)
                G = Lines_props[line][:g_total]
                B = Lines_props[line][:b_total]

                # EQ (36) - (37) #
                # active flow
                @constraint(model,
                    active_flow[t, line] == options.Vnom * (G * (nodal_voltage[t, NODE_FROM[idx]] - nodal_voltage[t, NODE_TO[idx]]) + 
                                            options.Vnom * B * (nodal_angle[t, NODE_TO[idx]] - nodal_angle[t, NODE_FROM[idx]]))
                )

                # reactive flow
                @constraint(model,
                    reactive_flow[t, line] == options.Vnom * (B * (nodal_voltage[t, NODE_TO[idx]] - nodal_voltage[t, NODE_FROM[idx]]) + 
                                            options.Vnom * G * (nodal_angle[t, NODE_TO[idx]] - nodal_angle[t, NODE_FROM[idx]]))
                )

                # voltage angle difference limits (in degrees)
                # assumed to be limited by 30 deg
                @constraint(model,
                    nodal_angle[t, NODE_FROM[idx]] - nodal_angle[t, NODE_TO[idx]] ≤ 30
                )

                @constraint(model,
                    nodal_angle[t, NODE_FROM[idx]] - nodal_angle[t, NODE_TO[idx]] ≥ -30
                )

                # EQ (42) #
                # linearised thermal constraints
                @constraint(model,
                    active_flow[t, line] + reactive_flow[t, line] ≤ sqrt(2) * Lines_props[line][:s_max]
                )

                @constraint(model,
                    active_flow[t, line] - reactive_flow[t, line] ≤ sqrt(2) * Lines_props[line][:s_max]
                )

                @constraint(model,
                    -active_flow[t, line] + reactive_flow[t, line] ≤ sqrt(2) * Lines_props[line][:s_max]
                )

                @constraint(model,
                    -active_flow[t, line] - reactive_flow[t, line] ≤ sqrt(2) * Lines_props[line][:s_max]
                )
            end
        end

    elseif options.power_flow == :dc

    #TODO: add DC power flow equations
    
    end

end     # end power_flow_Constraints
