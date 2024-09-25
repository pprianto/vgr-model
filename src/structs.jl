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
    El_Heat_Tax     ::Float64 = 60.0
    Vnom            ::Float64 = 130.0   # system voltage level
    FlexLim         ::Symbol = :yes
    EV              ::Symbol = :yes
    ### others coming up
end


Base.@kwdef mutable struct EVOptions
    BatCap          ::Float64 = 60      # assumed battery capacity, default 60 kWh (15, 30, 60, 85)
    CP              ::Float64 = 6.9     # assumed charger power, default 6.9 kW (3.7, 6.9, 11.0, 22.0)
    Charging_infra  ::Symbol = :h1      # assumed passenger car charging infrastructure (home, h1, h3, h6, ers)
    V2G             ::Symbol = :no      # charging capability :yes or :no
    Direct          ::Symbol = :no      # charging capability :yes or :no
    Optimal         ::Symbol = :no      # charging capability :yes or :no
    Qty             ::Symbol = :high    # scenario of number of EVs :high or :low
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


mutable struct ModelVariables # Variables
    total_cost              ::VariableRef
    capex                   ::VariableRef
    fix_om                  ::VariableRef
    fuel_cost               ::VariableRef
    var_om                  ::VariableRef
    start_part_costs        ::VariableRef
    exp_imp_costs           ::VariableRef
    tax_cost                ::VariableRef
    existing_generation     ::JuMP.Containers.DenseAxisArray
    generation_investment   ::JuMP.Containers.DenseAxisArray
    storage_investment      ::JuMP.Containers.DenseAxisArray
    active_generation       ::JuMP.Containers.DenseAxisArray
    reactive_generation     ::JuMP.Containers.DenseAxisArray
    generation_spin         ::Union{JuMP.Containers.DenseAxisArray, Nothing}
    generation_on           ::Union{JuMP.Containers.DenseAxisArray, Nothing}
    gen_startup_cost        ::Union{JuMP.Containers.DenseAxisArray, Nothing}
    gen_partload_cost       ::Union{JuMP.Containers.DenseAxisArray, Nothing}
    gen_startup_CO2         ::Union{JuMP.Containers.DenseAxisArray, Nothing}
    gen_partload_CO2        ::Union{JuMP.Containers.DenseAxisArray, Nothing}
    storage_charge          ::JuMP.Containers.DenseAxisArray
    storage_discharge       ::JuMP.Containers.DenseAxisArray
    storage_level           ::JuMP.Containers.DenseAxisArray
    nodal_voltage           ::JuMP.Containers.DenseAxisArray
    nodal_angle             ::JuMP.Containers.DenseAxisArray
    import_export           ::JuMP.Containers.DenseAxisArray
    export_to               ::JuMP.Containers.DenseAxisArray
    import_from             ::JuMP.Containers.DenseAxisArray
    active_flow             ::JuMP.Containers.DenseAxisArray
    reactive_flow           ::JuMP.Containers.DenseAxisArray
    pev_charging_slow       ::Union{JuMP.Containers.DenseAxisArray, Nothing}
    pev_discharge_net       ::Union{JuMP.Containers.DenseAxisArray, Nothing}
    pev_storage             ::Union{JuMP.Containers.DenseAxisArray, Nothing}
    pev_need                ::Union{JuMP.Containers.DenseAxisArray, Nothing}
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

