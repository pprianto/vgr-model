# function ev_agg(
#     model::Model,
#     sets::ModelSets,
#     params::ModelParameters,
#     vars::ModelVariables,
#     grid_infra::GridInfrastructures,
#     profiles::Profiles
# )

#=------------------------------------------------------------------------------
----- EV ADD on to MULTINODE
----- Year: 2021
----- by: Maria Taljeg�rd
----- last modified: 210421
-------------------Scenarier i ELIN----------------------------------------------------------------------------------------------------------------------------------------------------------------
The main unit is [MWh]
But there could be numerical difficulties with some model runs. Then a suggestions is to change the unit of the VPCB_discharge to MWh or VPCB_need to kWh.
the model can be run for all European countries. Number of electric vehicles is country dependent. Driving pattern and distance is based on the data from Sweden.
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Passenger cars only. direct = charging directly when being connected to the grid after a trip, 
optimal = optimised the charging while prioritising trips, V2G= same as optimal but includes also discharging to the grid
H2_EUroadmap then a share of the transport is running on hydrogen based on numbers from the EU hydrogen roadmap and a share on electricity
one can select multiple of the options direct/optimal/V2G and let a share of the fleet use one of the strategies
------------------------------------------------------------------------------=#
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


# other road vehicle types than passenger cars, "LT" = light trucks, "HT" = heavy trucks and "Bus"=bus
TRUCKS_BUSES = [
    "LT",
    "HT",
    "BUS"
]

EV_eff = 0.95

# Summarised charging infrastructure in Dict for easier calling
EV_infrastructure = Dict(
    15 => Dict("EVBat_Cap" => 15, :home => 0.71, :h6 => 0.74, :h3 => 0.81, :h1 => 0.87, :ers => 1.0),
    30 => Dict("EVBat_Cap" => 30, :home => 0.83, :h6 => 0.91, :h3 => 0.92, :h1 => 0.95, :ers => 1.0),
    60 => Dict("EVBat_Cap" => 60, :home => 0.92, :h6 => 0.96, :h3 => 0.96, :h1 => 0.97, :ers => 1.0),
    85 => Dict("EVBat_Cap" => 85, :home => 0.94, :h6 => 0.97, :h3 => 0.99, :h1 => 0.99, :ers => 1.0)
)

# Fleet Battery Capacity (in kWh)
EVBat_Cap = EV_infrastructure[options.BatCap]["EVBat_Cap"]

# Charging Power (kW)
if EV.CP == 3.7
    CP_slow = 3.7
elseif EV.CP == 6.9
    CP_slow = 6.9
elseif EV.CP == 11.0
    CP_slow = 11.0
elseif EV.CP == 22.0
    CP_slow = 22.0
end

# Share of V2G and optimum of the vehicle fleet
V2G_share = 0.3 # as default
opt_share = 0.3 # as default

# fuel (el or H2) consumption in kWh per km including losses from the engine in the vehicle
FC_el = 0.16
FC_H2 = 0.59*0.92*0.69

# estimated number of kilometers driven per year on el and non-el (in km)
km_yr = 13000

# Share of the kilometer for the vehicle fleet that can be driven on electricity
# depending mainly on charging infrastructure and battery size

share_el = EV_infrastructure[options.BatCap][options.Charging_infra]

@variables model begin
    pev_charging_slow[t ∈ PERIODS, i ∈ NODES] ≥ 0   # charging of the vehicle battery [MWh per hour]
    pev_discharge_net[t ∈ PERIODS, i ∈ NODES] ≥ 0   # Discharging of the vehicle battery back to the electricity grid [kWh per hour]
    pev_storage[t ∈ PERIODS, i ∈ NODES] ≥ 0         # Storage level of the vehicle battery [MWh per hour]
    pev_need[t ∈ PERIODS, i ∈ NODES] ≥ 0            # vehicle kilometers not met by charging [MWh per hour]
end

# Current initial storage level assumed to be 0
initEVSto = zeros(length(NODES))#, length(STO_TECHS))#, length(PERIODS))
Initial_EV_Storage = AxisArrays.AxisArray(
                                initEVSto, 
                                AxisArrays.Axis{:node_id}(NODES), 
                                # AxisArrays.Axis{:Tech}(STO_TECHS)
)#, Axis{:time}(PERIODS)) 

@constraint(model, EV_storage_level[t ∈ PERIODS, i ∈ NODES],
    pev_storage[t, i] ==
    (t > 1 ? pev_storage[t-1, i] : Initial_EV_Storage[i]) + 
    pev_charging_slow[t, i] * EV_eff +
    pev_need[t, i] -
    pev_discharge_net[t, i] #-
    #TODO: demand as parameter
)


# end