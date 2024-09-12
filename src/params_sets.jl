function make_sets(
    price::Prices,
    grid_infra::GridInfrastructures,
    tech_props::TechProps,
    demand::Demands,
    # run::Symbol
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
        Eldemand_data = demand.el[1:24, :]
        Heatdemand_data = demand.heat[1:24, :]
        H2demand_data = demand.h2[1:24, :]

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

    Gentech_data = Dict(tech_props.gen[!, :Tech] .=> eachrow(tech_props.gen[!, Not(:Tech)]))
    Stotech_data = Dict(tech_props.sto[!, :Tech] .=> eachrow(tech_props.sto[!, Not(:Tech)]))

    #=------------------------------------------------------------------------------
    POWER FACTOR
    ------------------------------------------------------------------------------=#
    # Power factor assumptions
    Gen_cos_ϕ = [0.8, 1]                            # assumed power factor operation bounds of generation
    Gen_sin_ϕ = sqrt.(1 .- Gen_cos_ϕ.^2)
    Demand_cos_ϕ = 0.95                             # assumed load power factor
    Demand_sin_ϕ = sqrt.(1 .- Demand_cos_ϕ.^2)
    
    # Assuming the reactive demand corresponds to 0.95 cos phi
    Reactive_Demand = Eldemand_data .* Demand_sin_ϕ ./ Demand_cos_ϕ 

    #=------------------------------------------------------------------------------
    POWER LINES PROPERTIES
    ------------------------------------------------------------------------------=#
    # Extract power lines sets
    # and properties
    lines_df, Ybus, G_Ybus, B_Ybus = lines_props(grid_infra.lines)
    LINES_SETS, Lines_props = arcs_prep(lines_df)

    @unpack LINES, 
            NODE_FROM, 
            NODE_TO = LINES_SETS

    #=------------------------------------------------------------------------------
    MODEL SETS
    ------------------------------------------------------------------------------=#

    #=------------------------------------------------------------------------------
    MAIN SETS
    ------------------------------------------------------------------------------=#  
    NODES = grid_infra.subs[!, :node_id]        # node set
    GEN_TECHS = tech_props.gen[!, :Tech]        # generation tech set
    STO_TECHS = tech_props.sto[!, :Tech]        # storage tech set

    if options.run == :full
        PERIODS = demand.el[!, :hour]               # time period set (hourly), full run
    elseif options.run == :trial || :test
        PERIODS = demand.el[1:24, :hour]             # time period set (hourly), trial runs
    else
        @error "No run type named $run."
    end
    #=------------------------------------------------------------------------------
    SUBSETS
    ------------------------------------------------------------------------------=#     
    # Coastal municipalities
    # for offshore wind eligibility
    # Strömstad
    # Tanum
    # Lysekil
    # Orust
    # Stenungsund
    # Kungälv
    # Göteborg
    COAST_NODES = [
        "GBG1",     # 1
        "GBG2",
        "GBG3",
        "GBG4",
        "GBG5",
        "GBG6",
        "GBG7",
        "GBG8",
        "GBG9",
        "GBG10",    # 10
        "GBG11",
        "GBG12",
        "KUN1",
        "KUN2",
        "KUN3",
        "KUN4",
        "STE1",
        "ORU1",
        "ORU2",
        "LYS1",     # 20
        "LYS2",
        "TAN1",
        "TAN2",
        "TAN3",
        "TAN4",
        "STR1",
        "STR2",
        "STR3",     # 28
    ]

    # Pit Thermal Storage is not feasible in Gothenburg
    GBG = [
        "GBG1",     # 1
        "GBG2",
        "GBG3",
        "GBG4",
        "GBG5",
        "GBG6",
        "GBG7",
        "GBG8",
        "GBG9",
        "GBG10",    # 10
        "GBG11",
        "GBG12",
    ]

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
            "COCHP",    # coal chp
            "WCHP",     # waste chp
            "WCCHP",    # wood chips chp
            "WPCHP",    # wood pellets chp
            # "SBCHP"     # straw biomass chp
    ]

    FC = [              # fuel cells
            "SOFC",     # solid oxide fuel cells
            # "PFC"       # pem fuel cell
    ]

    WIND = [
            "WON",      # onshore
            "WOFF"      # offshore
    ]

    PV = [
            "PVROOF",   # rooftop pv residential
            "PVUTIL",   # fixed axis utility PV
            "PVTRACK"   # single axis tracking utility PV
    ]

    HP = [
            "HPAIR",    # heat pumps air source
            "HPEX",      # heat pumps excess heat
            "HPSW"      # seawater heat pumps
    ]

    BOILER = [
                "EB",   # electric boilers
                "GB",    # natural gas boilers
                "HBW"   # biomass heat only boiler
    ]
    
    EC = [
            "AEC",      # alkali electrolyser
            "PEMEC"     # pemec electrolyser
    ]

    # thermal techs where flex lim applies
    FLEX_TH = [
                "COCHP",
                "GTSC",
                "CCGT",
                "GEBG",
                "WCHP",
                "WCCHP",
                "WPCHP",
                # "SBCHP"
    ]

    # techs with 1 PERIODS or less start up time
    THERMAL_1H = [
                    "GTSC",
                    "GEBG"
    ]

    # techs with 2 PERIODS start up time
    THERMAL_2H = [
                    "CCGT",
                    "WCHP"
    ]

    # # techs with 8 PERIODS start up time
    # THERMAL_8H = [
    #                 "SBCHP"
    # ]

    # techs with 12 PERIODS start up time
    THERMAL_12H = [
                    "COCHP",
                    "WCCHP",
                    "WPCHP"
    ]

    #=------------------------------------------------------------------------------
    Return
    ------------------------------------------------------------------------------=#

    sets =  ModelSets( 
            NODES, 
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
            THERMAL_12H
    )

    params = ModelParameters( 
             SE3_price,
             NO1_price,
             DK1_price, 
             Gentech_data, 
             Stotech_data, 
             Eldemand_data,
             Reactive_Demand, 
             Heatdemand_data, 
             H2demand_data, 
            #  Discount_rate,
             Gen_cos_ϕ, 
             Gen_sin_ϕ, 
             Demand_cos_ϕ, 
             Demand_sin_ϕ, 
             Lines_props,
            #  Vnom 
    )

    return  sets, params

end
