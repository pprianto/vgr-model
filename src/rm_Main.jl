#=------------------------------------------------------------------------------
------------------------------- REGION MODEL -----------------------------------

Main script to run the model:
include the following steps:
    1. reading input data
    2. initiate model
    3. run model
    4. ...

Return:
    - Dataframe of results for ...
------------------------------------------------------------------------------=#

using DataFrames, AxisArrays, CSV, XLSX, UnPack
using JuMP, HiGHS, Ipopt, GLPK, Gurobi
import SparseArrays, LinearAlgebra

include("rm_Utilities.jl")
include("rm_Params_Sets.jl")
include("rm_Vars_Consts.jl")
include("rm_Model.jl")

# Set directory
current_dir = pwd()
mkpath("modelinput")        # input folder
mkpath("results")           # results folder
input_dir = joinpath(current_dir, "modelinput")
results_dir = joinpath(current_dir, "results")

#=---------------------------------------------
INPUT DATA
---------------------------------------------=#

# electricity price
elpris_df = read_file(joinpath(input_dir, "elpris_.xlsx"))

# electrical infrastructure
substations_df = read_file(joinpath(input_dir, "subs_final.csv"))
lines_df = read_file(joinpath(input_dir, "lines_final.csv"))
pp_df = read_file(joinpath(input_dir, "pp_final.csv"))
pp_df = rename_pp(pp_df)                                                                    # rename the tech according to the index sets

# technology properties
# assume 2050
gen_tech_df = read_file(joinpath(input_dir, "gen_tech_2050_.xlsx"))
sto_tech_df = read_file(joinpath(input_dir, "sto_tech_2050_.xlsx"))

# demand data
# currently in hourly period
el_demand_df = read_file(joinpath(input_dir, "el_nodal_demand.csv"))
heat_demand_df = read_file(joinpath(input_dir, "heat_nodal_demand.csv"))
h2_demand_df = read_file(joinpath(input_dir, "h2_nodal_demand.csv"))

# RE profile
PV_fix_profile = read_file(joinpath(input_dir, "nodal_profile_pv_fixed_2019.csv"))          # fixed axis pv
PV_opt_profile = read_file(joinpath(input_dir, "nodal_profile_pv_double_axis_2019.csv"))    # opt tracking pv
WT_on_profile = read_file(joinpath(input_dir, "nodal_profile_onshore_wind_2019.csv"))       # onshore wind
WT_off_profile = read_file(joinpath(input_dir, "nodal_profile_offshore_wind_2019.csv"))     # offshore wind

# collections of input data in NamedTuple
price = (; el=elpris_df)
grid_infra = (; subs=substations_df, lines=lines_df, pp=pp_df)
tech_props = (; gen=gen_tech_df, sto=sto_tech_df)
demand = (; el=el_demand_df, heat=heat_demand_df, h2=h2_demand_df)
profiles = (; PVfix=PV_fix_profile, PVopt=PV_opt_profile, WTon=WT_on_profile, WToff=WT_off_profile)

#=---------------------------------------------
INITIATE MODEL
---------------------------------------------=#
model = Model(Gurobi.Optimizer)
# set_attribute(model, "Presolve", 0)
# set_optimizer_attribute(model, "NumericFocus", 2)
# set_optimizer_attribute(model, "BarHomogeneous", 1)     # due to rhs? TODO: improve the area for invesment constraint
set_optimizer_attribute(model, "Threads", Threads.nthreads())
log_file = joinpath(results_dir, "model_log.txt")
set_optimizer_attribute(model, "LogFile", log_file)

@time model, results, times = run_Model(
    model,
    price,
    grid_infra,
    tech_props,
    demand,
    profiles  
);

