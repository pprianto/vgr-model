function make_sets(
    price::Prices,
    grid_infra::GridInfrastructures,
    tech_props::TechProps,
    demand::Demands,
)

#=------------------------------------------------------------------------------
----------------------------- SETS & PARAMETERS --------------------------------

I. Define Sets of the model
current sets:
1. nodes (i)
2. generation techs (x)
3. storage techs (s)
4. times (t) - hourly
5. ...

II. Define Parameters of the model
current parameters:
1. cost parameters for techs
2. electricity price for import
3. demand (el and heat)
4. tech props
5. ...


Return
1. Sets NamedTuple

2. Parameters NamedTuple

------------------------------------------------------------------------------=#
    #=------------------------------------------------------------------------------
    MODEL PARAMETERS / INPUT DATA
    ------------------------------------------------------------------------------=#

    #=------------------------------------------------------------------------------
    DEMAND
    ------------------------------------------------------------------------------=#
    if options.run == :full
        # demand data in each nodes
        # full run
        Eldemand_data = demand.el
        Heatdemand_data = demand.heat
        H2demand_data = demand.h2

    elseif options.run == :trial || :test
        # trial runs
        # adjust hour accordingly
        Eldemand_data = demand.el[1:168, :]
        Heatdemand_data = demand.heat[1:168, :]
        H2demand_data = demand.h2[1:168, :]

    else
        @error "No run type named $run."
    
    end
    #=------------------------------------------------------------------------------
    PRICES
    ------------------------------------------------------------------------------=#
    # electricity prices in €/MWh 
    SE3_price = price.SE3
    NO1_price = price.NO1
    DK1_price = price.DK1         

    # heat and hydrogen market price?
    # fuel price? 

    #=------------------------------------------------------------------------------
    GENERATION AND STORAGE TECHNOLOGY PROPERTIES
    ------------------------------------------------------------------------------=#
    # Dict of tech properties, keys are tech names, values are properties in df format
    # based on Danish Energy Agency catalogues
    # https://ens.dk/en/our-services/technology-catalogues/technology-data-generation-electricity-and-district-heating
    # https://ens.dk/en/our-services/technology-catalogues/technology-data-energy-storage
    # https://ens.dk/en/our-services/technology-catalogues/technology-data-renewable-fuels
    # https://ens.dk/en/our-services/technology-catalogues/technology-data-transport-energy

    Gentech_data = Dict(Symbol.(tech_props.gen[!, :Tech]) .=> NamedTuple.(eachrow(tech_props.gen[!, Not(:Tech)])))
    Stotech_data = Dict(Symbol.(tech_props.sto[!, :Tech]) .=> NamedTuple.(eachrow(tech_props.sto[!, Not(:Tech)])))

    #=------------------------------------------------------------------------------
    POWER FACTOR
    ------------------------------------------------------------------------------=#
    # Power factor assumptions
    # Gen_cos_ϕ = [0.8, 1]                            # assumed power factor operation bounds of generation
    # Gen_sin_ϕ = sqrt.(1 .- Gen_cos_ϕ.^2)
    Demand_cos_ϕ = 0.95                             # assumed load power factor
    Demand_sin_ϕ = sqrt.(1 .- Demand_cos_ϕ.^2)
    
    # Assuming the reactive demand corresponds to 0.95 cos phi
    Reactive_Demand = Eldemand_data .* Demand_sin_ϕ ./ Demand_cos_ϕ 

    #=------------------------------------------------------------------------------
    POWER LINES PROPERTIES
    ------------------------------------------------------------------------------=#
    # Extract power lines sets
    # and properties
    prepped_lines = lines_props(grid_infra.lines)
    LINES_SETS, Lines_props = lines_prep(prepped_lines.lines_df)

    (;  LINES, 
        NODE_FROM, 
        NODE_TO
    ) = LINES_SETS

    #=------------------------------------------------------------------------------
    MODEL SETS
    ------------------------------------------------------------------------------=#

    #=------------------------------------------------------------------------------
    MAIN SETS
    ------------------------------------------------------------------------------=#  
    NODES = Symbol.(grid_infra.subs[!, :node_id])        # node set
    GEN_TECHS = Symbol.(tech_props.gen[!, :Tech])        # generation tech set
    STO_TECHS = Symbol.(tech_props.sto[!, :Tech])        # storage tech set

    if options.run == :full
        PERIODS = demand.el[!, :hour]               # time period set (hourly), full run
    elseif options.run == :trial || :test
        PERIODS = demand.el[1:168, :hour]             # time period set (hourly), trial runs
    else
        @error "No run type named $run."
    end
    #=------------------------------------------------------------------------------
    SUBSETS
    ------------------------------------------------------------------------------=#     
    # Coastal municipalities
    # for offshore wind and sea water hp eligibility
    # Strömstad
    # Tanum
    # Lysekil
    # Orust
    # Stenungsund
    # Kungälv
    # Göteborg

    coast_municipalities = ["GBG"; "KUN"; "STE"; "ORU"; "LYS"; "TAN"; "STR"]
    coast_subs_df = filter(row -> any(startswith(row.node_id, prefix) for prefix in coast_municipalities), grid_infra.subs)
    
    COAST_NODES = Symbol.(coast_subs_df.node_id |> Vector)

    # Pit Thermal Storage is not feasible in Gothenburg

    gothenburg = ["GBG"]
    gbg_subs_df = filter(row -> any(startswith(row.node_id, prefix) for prefix in gothenburg), grid_infra.subs)
    
    GBG = Symbol.(gbg_subs_df.node_id |> Vector)

    #=------------------------------------------------------------------------------
    GENERATION TECHNOLOGY SUBSETS
    ------------------------------------------------------------------------------=# 
    # carrier generation technologies subset
    EL_GEN, HEAT_GEN, H2_GEN = carrier_subsets(tech_props.gen) 

    # carrier storage technologies subset
    EL_STO, HEAT_STO, H2_STO = carrier_subsets(tech_props.sto)    

    #=------------------------------------------------------------------------------
    TRANSMISSION NODES SUBSETS
    ------------------------------------------------------------------------------=# 
    trans_node = filter(row -> row.import_trans == true, grid_infra.subs)    
    SE3_TRANS_NODES = filter(row -> !(row.node_id in ["MOL1", "DAL3"]), trans_node).node_id
    NO1_TRANS_NODES = filter(row -> row.node_id == "DAL3", trans_node).node_id
    DK1_TRANS_NODES = filter(row -> row.node_id == "MOL1", trans_node).node_id
    TRANSMISSION_NODES = trans_node.node_id

    #=------------------------------------------------------------------------------
    SPECIFIC TECHNOLOGY SUBSETS
    ------------------------------------------------------------------------------=# 
    # manually defined from excel file input
    CHP = [             # CHP
            :COCHP,    # coal chp
            :WCHP,     # waste chp
            :WCCHP,    # wood chips chp
            :WPCHP,    # wood pellets chp
            # :SBCHP     # straw biomass chp
    ]

    FC = [              # fuel cells
            :SOFC,     # solid oxide fuel cells
            # :PFC:       # pem fuel cell
    ]

    WIND = [
            :WON,      # onshore
            :WOFF      # offshore
    ]

    PV = [
            :PVROOF,   # rooftop pv residential
            :PVUTIL,   # fixed axis utility PV
            :PVTRACK   # single axis tracking utility PV
    ]

    HP = [
            :HPAIR,    # heat pumps air source
            :HPEX,      # heat pumps excess heat
            :HPSW      # seawater heat pumps
    ]

    BOILER = [
                :EB,   # electric boilers
                :GB,    # natural gas boilers
                :HBW   # biomass heat only boiler
    ]
    
    EC = [
            :AEC,      # alkali electrolyser
            :PEMEC     # pemec electrolyser
    ]

    # thermal techs where flex lim applies
    FLEX_TH = [
                :COCHP,
                :GTSC,
                :CCGT,
                :GEBG,
                :WCHP,
                :WCCHP,
                :WPCHP,
                # :SBCHP
    ]

    # techs with 1 PERIODS or less start up time
    THERMAL_1H = [
                    :GTSC,
                    :GEBG
    ]

    # techs with 2 PERIODS start up time
    THERMAL_2H = [
                    :CCGT,
                    :WCHP
    ]

    # # techs with 8 PERIODS start up time
    # THERMAL_8H = [
    #                 :SBCHP:
    # ]

    # techs with 12 PERIODS start up time
    THERMAL_12H = [
                    :COCHP,
                    :WCCHP,
                    :WPCHP
    ]

    # Batteries separated into energy and capacity components
    BAT_CAP = [
                :LI_CAP; 
                :VR_CAP
    ]
    
    BAT_EN = [
                :LI_EN; 
                :VR_EN
    ]

    # subsets of storages, excluding capacity components of batteries
    STO_EN = setdiff(STO_TECHS, BAT_CAP)

    #=------------------------------------------------------------------------------
    Return
    ------------------------------------------------------------------------------=#

    sets =  (;  NODES = Symbol.(NODES), 
                TRANSMISSION_NODES = Symbol.(TRANSMISSION_NODES),
                SE3_TRANS_NODES = Symbol.(SE3_TRANS_NODES),
                NO1_TRANS_NODES = Symbol.(NO1_TRANS_NODES),
                DK1_TRANS_NODES = Symbol.(DK1_TRANS_NODES),
                COAST_NODES = Symbol.(COAST_NODES),
                GBG = Symbol.(GBG),
                GEN_TECHS = Symbol.(GEN_TECHS), 
                EL_GEN = Symbol.(EL_GEN), 
                HEAT_GEN = Symbol.(HEAT_GEN), 
                H2_GEN = Symbol.(H2_GEN), 
                STO_TECHS = Symbol.(STO_TECHS), 
                EL_STO = Symbol.(EL_STO), 
                HEAT_STO = Symbol.(HEAT_STO), 
                H2_STO = Symbol.(H2_STO), 
                PERIODS = Array(PERIODS), 
                LINES = Symbol.(LINES), 
                NODE_FROM = Symbol.(NODE_FROM), 
                NODE_TO = Symbol.(NODE_TO), 
                CHP = Symbol.(CHP),
                FC = Symbol.(FC),
                WIND = Symbol.(WIND),
                PV = Symbol.(PV),
                HP = Symbol.(HP),
                BOILER = Symbol.(BOILER),
                EC = Symbol.(EC),
                FLEX_TH = Symbol.(FLEX_TH),
                THERMAL_1H = Symbol.(THERMAL_1H),
                THERMAL_2H = Symbol.(THERMAL_2H),
                # THERMAL_8H = Symbol.(THERMAL_8H),
                THERMAL_12H = Symbol.(THERMAL_12H),
                BAT_CAP = Symbol.(BAT_CAP),
                BAT_EN = Symbol.(BAT_EN),
                STO_EN = Symbol.(STO_EN)
    )

    params = (; SE3_price,
                NO1_price,
                DK1_price, 
                Gentech_data, 
                Stotech_data, 
                Eldemand_data,
                Reactive_Demand, 
                Heatdemand_data, 
                H2demand_data, 
                # Discount_rate,
                # Gen_cos_ϕ, 
                # Gen_sin_ϕ, 
                Demand_cos_ϕ, 
                Demand_sin_ϕ, 
                Lines_props,
                # Vnom 
    )

    return  sets, params

end


function make_variables(
    model::Model,
    # sets::ModelSets,
    # params::ModelParameters
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

    (;  NODES, 
        TRANSMISSION_NODES,
        GEN_TECHS, 
        EL_GEN, 
        # HEAT_GEN, 
        # H2_GEN, 
        STO_TECHS, 
        PERIODS, 
        LINES, 
        FLEX_TH,
        STO_EN,
        HP,
        EC,
        FC,
    ) = sets
        
    (; Lines_props) = params
            # Vnom = params

    #=------------------------------------------------------------------------------
    MODEL VARIABLES
    ------------------------------------------------------------------------------=#

    # Cost-related variables
    # all free variables
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
        active_generation[t ∈ PERIODS, i ∈ NODES, x ∈ GEN_TECHS]  ≥ 0
        reactive_generation[t ∈ PERIODS, i ∈ NODES, x ∈ EL_GEN]
    end

    # temporary variables for conversion technologies
    @variables model begin
        HP_as_demand[t ∈ PERIODS, i ∈ NODES, x ∈ HP] ≥ 0    # Heat Pumps
        EB_as_demand[t ∈ PERIODS, i ∈ NODES, x ∈ [:EB]] ≥ 0       # Electric Boiler
        EC_as_demand[t ∈ PERIODS, i ∈ NODES, x ∈ EC] ≥ 0    # Electrolyser
        FC_as_demand[t ∈ PERIODS, i ∈ NODES, x ∈ FC] ≥ 0    # Fuel Cell
    end

    # Storage-related variables (MWh)
    @variables model begin
        storage_charge[t ∈ PERIODS, i ∈ NODES, s ∈ STO_EN]     ≥ 0
        storage_discharge[t ∈ PERIODS, i ∈ NODES, s ∈ STO_EN]  ≥ 0
        storage_level[t ∈ PERIODS, i ∈ NODES, s ∈ STO_EN]      ≥ 0
    end

    #=------------------------------------------------------------------------------
    TWO-VARIABLE APPROACH VARIABLES
    ------------------------------------------------------------------------------=#
    # applies to thermal power plants
    # related to start up / ramping limits
    # based on Lisa's thesis
    # only applicable if Flexlim option is enabled
    if options.FlexLim == :yes
        @variables model begin
            generation_spin[t ∈ PERIODS, i ∈ NODES, x ∈ FLEX_TH]      ≥ 0
            generation_on[t ∈ PERIODS, i ∈ NODES, x ∈ FLEX_TH]        ≥ 0
            gen_startup_cost[t ∈ PERIODS, i ∈ NODES, x ∈ FLEX_TH]     ≥ 0
            gen_partload_cost[t ∈ PERIODS, i ∈ NODES, x ∈ FLEX_TH]    ≥ 0
            gen_startup_CO2[t ∈ PERIODS, i ∈ NODES, x ∈ FLEX_TH]      ≥ 0
            gen_partload_CO2[t ∈ PERIODS, i ∈ NODES, x ∈ FLEX_TH]     ≥ 0
        end
    end

    #=------------------------------------------------------------------------------
    EXPORT - IMPORT VARIABLES
    EQ (38)
    ------------------------------------------------------------------------------=#
    # EQ (33) - (35) #
    # export and import to/from transmission system
    # from ACCEL report, the subscription of Vattenfall Eldistribution - SvK is ~3200MW ≈ 400MW per nodes
    # perhaps trading regulations or annual limit?
    @variables model begin
            import_export[t ∈ PERIODS, i ∈ TRANSMISSION_NODES]    ≥ 0
        0 ≤ export_to[t ∈ PERIODS, i ∈ TRANSMISSION_NODES]        ≤ 400
        0 ≤ import_from[t ∈ PERIODS, i ∈ TRANSMISSION_NODES]      ≤ 400
    end

    # import-export constraints
    # note: import_export is used to simplify nodal balance and objective function
    # when import_export ≥ 0, then the nodes are importing and imposed to trade cost (elprice)
    # otherwise, then the nodes are exporting and gain from trade cost (elprice)
    @constraint(model, 
                Trade[t ∈ PERIODS, i ∈ TRANSMISSION_NODES],
                import_export[t, i] == import_from[t, i] - export_to[t, i]
    )

    #=------------------------------------------------------------------------------
    VOLTAGE VARIABLES
    EQ (38)
    ------------------------------------------------------------------------------=#
    # voltage nominal and angles (in kV for voltage magnitude, degree for angle)
    # first nodes is excluded since it is considered slack bus
    @variables model begin
        0.9 * options.Vnom ≤ nodal_voltage[t ∈ PERIODS, i ∈ NODES] ≤ 1.1 * options.Vnom 
        -360 ≤ nodal_angle[t ∈ PERIODS, i ∈ NODES] ≤ 360
    end

    # Slack bus voltage and angle over time, assumes the first entry as slack

    for t ∈ PERIODS
        fix(nodal_voltage[t, NODES[1]], 
            options.Vnom,
            force=true
        )

        fix(nodal_angle[t, NODES[1]], 
            0.0,
            force=true
        )

    end

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
        -Lines_props[l][:s_max] ≤ active_flow[t ∈ PERIODS, l ∈ LINES]     ≤ Lines_props[l][:s_max]
        -Lines_props[l][:s_max] ≤ reactive_flow[t ∈ PERIODS, l ∈ LINES]   ≤ Lines_props[l][:s_max]
    end

    #=------------------------------------------------------------------------------
    EV-RELATED VARIABLES
    MARIA'S MODEL
    ------------------------------------------------------------------------------=#
    # only applicable if EV option is enabled
    if options.EV == :yes
        @variables model begin
            pev_charging_slow[t ∈ PERIODS, i ∈ NODES] ≥ 0   # charging of the vehicle battery [MWh per hour]
            pev_discharge_net[t ∈ PERIODS, i ∈ NODES] ≥ 0   # Discharging of the vehicle battery back to the electricity grid [kWh per hour]
            pev_storage[t ∈ PERIODS, i ∈ NODES] ≥ 0         # Storage level of the vehicle battery [MWh per hour]
            pev_need[t ∈ PERIODS, i ∈ NODES] ≥ 0            # vehicle kilometers not met by charging [MWh per hour]
        end
    end

    #=------------------------------------------------------------------------------
    Return
    ------------------------------------------------------------------------------=#

    if options.FlexLim == :yes && options.EV == :yes

        return  (;  total_cost,
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
                    reactive_flow,
                    pev_charging_slow,       
                    pev_discharge_net,       
                    pev_storage,             
                    pev_need,      
                    HP_as_demand,
                    EB_as_demand,
                    EC_as_demand,
                    FC_as_demand,          
        )

    elseif options.FlexLim == :yes && options.EV == :no

        return  (;  total_cost,
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
                    reactive_flow,
                    HP_as_demand,
                    EB_as_demand,
                    EC_as_demand,
                    FC_as_demand,   
        )

    else

        return  (;  total_cost,
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
                    storage_charge,
                    storage_discharge,
                    storage_level,
                    nodal_voltage,
                    nodal_angle,
                    import_export,
                    export_to,
                    import_from,
                    active_flow,
                    reactive_flow,
                    HP_as_demand,
                    EB_as_demand,
                    EC_as_demand,
                    FC_as_demand,   
        )

    end

end         # end make_Variables


function lines_prep(
    lines_df::DataFrame
)

    lines = lines_df[!, [:lines_id, :node_from, :node_to, :g_total, :b_total, :s_max]]
    lines[!, :arcs_fr] = [(row.lines_id, row.node_from, row.node_to) for row in eachrow(lines)]
    lines[!, :arcs_to] = [(row.lines_id, row.node_to, row.node_from) for row in eachrow(lines)]
    # lines[!, :arcs] = [lines[!, :arcs_fr]; lines[!, :arcs_to]]

    # Power system sets
    LINES = lines[!, :lines_id]         # power lines set
    NODE_FROM = lines[!, :node_from]    # set of node from of the line
    NODE_TO = lines[!, :node_to]        # set of node to of the line
    ARCS_FR = lines[!, :arcs_fr]        # combined set of lines - nodes
    ARCS_TO = lines[!, :arcs_to]

    Lines_props = Dict(Symbol.(lines[!, :lines_id]) .=> NamedTuple.(eachrow(lines[!, Not(:lines_id)])))

    lines_sets = (; 
                LINES,
                NODE_FROM, 
                NODE_TO, 
                ARCS_FR,
                ARCS_TO
            )

    return lines_sets, Lines_props

end

