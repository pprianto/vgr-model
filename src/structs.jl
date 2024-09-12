#=------------------------------------------------------------------------------
----------------------------------- STRUCTS ------------------------------------

Define structs to store model components parameters and results 

------------------------------------------------------------------------------=#

Base.@kwdef mutable struct ModelOptions
    # model options with pre-defined configs
    run             ::Symbol            # :full or :trial :test
    CO2_limit       ::Float64 = 0.0
    target_year     ::Int = 2050        # model target year
    profile_year    ::Int = 2019        # model profile
    power_flow      ::Symbol = :lac     # linearised AC (LAC) or DC
    Discount_rate   ::Float64 = 0.05    # assumed discount rate
    Vnom            ::Float64 = 130.0   # system voltage level

    ### others coming up
end


struct ModelSets # Sets
    NODES               ::Vector{String}
    TRANSMISSION_NODES  ::Vector{String}
    SE3_TRANS_NODES     ::Vector{String}
    NO1_TRANS_NODES     ::Vector{String}
    DK1_TRANS_NODES     ::Vector{String}
    COAST_NODES         ::Vector{String}
    GBG                 ::Vector{String}
    GEN_TECHS           ::Vector{String}
    EL_GEN              ::Vector{String}
    HEAT_GEN            ::Vector{String}
    H2_GEN              ::Vector{String}
    STO_TECHS           ::Vector{String}
    EL_STO              ::Vector{String}
    HEAT_STO            ::Vector{String}
    H2_STO              ::Vector{String}
    PERIODS             ::Vector{Int}
    LINES               ::Vector{String}
    NODE_FROM           ::Vector{String}
    NODE_TO             ::Vector{String}
    CHP                 ::Vector{String}
    FC                  ::Vector{String}
    WIND                ::Vector{String}
    PV                  ::Vector{String}
    HP                  ::Vector{String}
    BOILER              ::Vector{String}
    EC                  ::Vector{String}
    FLEX_TH             ::Vector{String}
    THERMAL_1H          ::Vector{String}
    THERMAL_2H          ::Vector{String}
    # THERMAL_8H          ::Vector{String}
    THERMAL_12H         ::Vector{String}
end


struct ModelParameters # Parameters
    SE3_price           ::Vector{Float64}
    NO1_price           ::Vector{Float64}
    DK1_price           ::Vector{Float64}
    Gentech_data        ::Dict{String, DataFrameRow}
    Stotech_data        ::Dict{String, DataFrameRow}
    Eldemand_data       ::DataFrame
    Reactive_Demand     ::DataFrame
    Heatdemand_data     ::DataFrame
    H2demand_data       ::DataFrame
    # Discount_rate       ::Float64
    Gen_cos_ϕ           ::Vector{Float64}
    Gen_sin_ϕ           ::Vector{Float64}
    Demand_cos_ϕ        ::Float64
    Demand_sin_ϕ        ::Float64
    Lines_props         ::Dict{String, DataFrameRow}
    # Vnom                ::Float64
end


struct ModelVariables # Variables
    total_cost          
    capex               
    fix_om              
    fuel_cost           
    var_om              
    start_part_costs    
    exp_imp_costs       
    tax_cost            
    existing_generation
    generation_investment 
    storage_investment 
    active_generation 
    reactive_generation 
    generation_spin
    generation_on
    gen_startup_cost
    gen_partload_cost
    gen_startup_CO2
    gen_partload_CO2            
    storage_charge 
    storage_discharge 
    storage_level 
    nodal_voltage 
    nodal_angle
    import_export 
    export_to 
    import_from 
    active_flow
    reactive_flow
end


struct ModelStruct
    model       ::Model             # model info
    sets        ::ModelSets         # model sets
    params      ::ModelParameters   # model parameters
    vars        ::ModelVariables    # model variables
    times       ::NamedTuple        # time to build and solve model
end


struct Prices
    SE3         ::Vector{Float64}
    NO1         ::Vector{Float64}
    DK1         ::Vector{Float64}
end


struct GridInfrastructures
    subs        ::DataFrame
    lines       ::DataFrame
    pp          ::DataFrame
end


struct TechProps
    gen         ::DataFrame
    sto         ::DataFrame
end


struct Demands
    el          ::DataFrame
    heat        ::DataFrame
    h2          ::DataFrame
end


struct Profiles
    PVfix       ::DataFrame
    PVopt       ::DataFrame
    WTon        ::DataFrame
    WToff       ::DataFrame
end

