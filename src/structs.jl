#=------------------------------------------------------------------------------
----------------------------------- STRUCTS ------------------------------------

Define structs to store model components parameters and results 
Also includes model options

------------------------------------------------------------------------------=#

Base.@kwdef mutable struct ModelOptions
    # model options with pre-defined configs
    run             ::Symbol            # :full or :trial :test
    CO2_budget      ::Float64 = 0.0     # in kgCO2
    target_year     ::Int = 2050        # model target year
    profile_year    ::Int = 2019        # model profile
    el_price_year   ::Int = 1991        # future el price model year (1991 or 1992)
    scenario        ::Symbol = :beta    # :now, :alpha, :beta, with beta as highest demand increase
    power_flow      ::Symbol = :lac     # linearised AC (LAC) or DC
    Discount_rate   ::Float64 = 0.05    # assumed discount rate
    El_Heat_Tax     ::Float64 = 60.0    # assumed tax to generate electricity, on top of el price (€/MWh)
    CO2_limit       ::Symbol = :fee     # constraint to limit CO2 emitting technology (:ton or :fee)
    Emission_fee    ::Float64 = 150.0   # assumed CO2 emission fee (€/tonCO2)
    Vnom            ::Float64 = 130.0   # system voltage level (kV)
    FlexLim         ::Symbol = :yes     # option if flexlim constraint is included
    EV              ::Symbol = :yes     # option if EV is included
    ### others coming up
end


Base.@kwdef mutable struct EVOptions
    BatCap          ::Float64 = 30      # assumed battery capacity, default 60 kWh (15, 30, 60, 85)
    CP              ::Float64 = 6.9     # assumed charger power, default 6.9 kW (3.7, 6.9, 11.0, 22.0)
    Charging_infra  ::Symbol = :h1      # assumed passenger car charging infrastructure (home, h1, h3, h6, ers)
    V2G             ::Symbol = :no      # charging capability :yes or :no
    Direct          ::Symbol = :yes      # charging capability :yes or :no
    Optimal         ::Symbol = :yes      # charging capability :yes or :no
    Qty             ::Symbol = :high    # scenario of number of EVs :high or :low
    ERS             ::Symbol = :no      # electrified road system
    Cars_H2         ::Symbol = :no
    Trucks_H2       ::Symbol = :no
    EV_number       ::Symbol = :high    # :high (60% 2030 100% 2050) or :low (30% 2030 60% 2050)
    ### others coming up
end


struct ModelSets # Sets
    NODES               ::Vector{Symbol}
    TRANSMISSION_NODES  ::Vector{Symbol}
    SE3_TRANS_NODES     ::Vector{Symbol}
    NO1_TRANS_NODES     ::Vector{Symbol}
    DK1_TRANS_NODES     ::Vector{Symbol}
    COAST_NODES         ::Vector{Symbol}
    GBG                 ::Vector{Symbol}
    GEN_TECHS           ::Vector{Symbol}
    EL_GEN              ::Vector{Symbol}
    HEAT_GEN            ::Vector{Symbol}
    H2_GEN              ::Vector{Symbol}
    STO_TECHS           ::Vector{Symbol}
    EL_STO              ::Vector{Symbol}
    HEAT_STO            ::Vector{Symbol}
    H2_STO              ::Vector{Symbol}
    PERIODS             ::Vector{Int}
    LINES               ::Vector{Symbol}
    NODE_FROM           ::Vector{Symbol}
    NODE_TO             ::Vector{Symbol}
    CHP                 ::Vector{Symbol}
    FC                  ::Vector{Symbol}
    WIND                ::Vector{Symbol}
    PV                  ::Vector{Symbol}
    HP                  ::Vector{Symbol}
    BOILER              ::Vector{Symbol}
    EC                  ::Vector{Symbol}
    FLEX_TH             ::Vector{Symbol}
    THERMAL_1H          ::Vector{Symbol}
    THERMAL_2H          ::Vector{Symbol}
    # THERMAL_8H          ::Vector{Symbol}
    THERMAL_12H         ::Vector{Symbol}
end


struct ModelParameters # Parameters
    SE3_price           ::Vector{Float64}
    NO1_price           ::Vector{Float64}
    DK1_price           ::Vector{Float64}
    Gentech_data        ::Dict{}
    Stotech_data        ::Dict{}
    Eldemand_data       ::DataFrame
    Reactive_Demand     ::DataFrame
    Heatdemand_data     ::DataFrame
    H2demand_data       ::DataFrame
    # Discount_rate       ::Float64
    # Gen_cos_ϕ           ::Vector{Float64}
    # Gen_sin_ϕ           ::Vector{Float64}
    Demand_cos_ϕ        ::Float64
    Demand_sin_ϕ        ::Float64
    Lines_props         ::Dict{}
    # Vnom                ::Float64
end


mutable struct ModelVariables # Variables
    # Union type for variables that can be excluded depending on model options
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

