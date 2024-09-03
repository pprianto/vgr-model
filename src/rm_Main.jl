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

using DataFrames, CSV, XLSX, UnPack, Printf
using JuMP, HiGHS, Ipopt, GLPK, Gurobi
import AxisArrays, SparseArrays, LinearAlgebra
import Plots, CairoMakie
CairoMakie.activate!()

include("rm_Utilities.jl")
include("rm_Params_Sets.jl")
include("rm_Vars_Consts.jl")
include("rm_Model.jl")

# Set directory
current_dir = pwd()
mkpath("modelinput")        # input folder, should already be created during retrive OSM process
mkpath("results")           # results folder
input_dir = joinpath(current_dir, "modelinput")
results_dir = joinpath(current_dir, "results")

#=---------------------------------------------
INPUT DATA
---------------------------------------------=#

# electricity price
elpris_df = read_file(joinpath(input_dir, "elpris.csv"))

# electrical infrastructure
substations_df = read_file(joinpath(input_dir, "subs_final.csv"))
lines_df = read_file(joinpath(input_dir, "lines_final.csv"))
pp_df = read_file(joinpath(input_dir, "pp_final.csv"))
pp_df = rename_pp(pp_df)                                                                    # rename the tech according to the index sets

# technology properties
# assume 2050
gen_tech_df = read_file(joinpath(input_dir, "gen_tech_2050.csv"))
sto_tech_df = read_file(joinpath(input_dir, "sto_tech_2050.csv"))

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
price = (; SE3=Array(elpris_df.SE3), NO1=Array(elpris_df.NO1), DK1=Array(elpris_df.DK1))
grid_infra = (; subs=substations_df, lines=lines_df, pp=pp_df)
tech_props = (; gen=gen_tech_df, sto=sto_tech_df)
demand = (; el=el_demand_df, heat=heat_demand_df, h2=h2_demand_df)
profiles = (; PVfix=PV_fix_profile, PVopt=PV_opt_profile, WTon=WT_on_profile, WToff=WT_off_profile)

#=---------------------------------------------
INITIATE MODEL
---------------------------------------------=#
# Gurobi parameters
optimizer = optimizer_with_attributes(
                Gurobi.Optimizer,
                "Threads" => Threads.nthreads(),
                "BarHomogeneous" => 1,              # 1: enabled
                "Crossover" => -1,                  # 0: disabled
                "Method" => -1,                     # -1: auto, 1: dual simplex
                "Presolve" => 2,                    # 2: aggressive
                # "NumericFocus" => 2,
)

m = Model(optimizer)
# set_silent(m)
# log_file = joinpath(results_dir, "model_log.txt")
# set_optimizer_attribute(m, "LogFile", log_file)

mutable struct Model_Struct
    model::Model        # model info
    sets::NamedTuple    # model sets
    params::NamedTuple  # model parameters
    vars::NamedTuple    # model variables
    times::NamedTuple   # time to build and solve model
end


@time model_done = run_Model(
    m,
    price,
    grid_infra,
    tech_props,
    demand,
    profiles  
);

(; model, sets, vars) = model_done


@time results = query_solutions(
    model,
    sets,
    vars
);

#=---------------------------------------------
POST - PROCESSING
---------------------------------------------=#

function vgr_investment(
    result_df
)

    if "NODE" in names(result_df)
        # get only VGR rows
        VGR_inv = result_df[result_df.NODE .== "VGR", 2:end]

        # remove anything that is zero (not investing)
        VGR_inv = VGR_inv[!, [col for col in names(VGR_inv) if VGR_inv[1, col] > 0]]

        return VGR_inv

    else
        println("column NODE is not found")

        return nothing
    
    end
end

vgr_gen_inv = vgr_investment(results.Generation_Investment)
vgr_sto_inv = vgr_investment(results.Storage_Investment)

function basic_barchart(
    inv_df
)

    x = names(inv_df)
    y = Array(inv_df[1,:])
    (ymin, ymax) = extrema(y)
    delta_y = 0.05*(ymax - ymin)
    cap_str = [(@sprintf("%.2f MWh", cap), 5, 45.0) for cap in y]    # the annotation, text size, text rotation

    Plots.bar(
        x,     
        y,                   
        legend = false,
        xticks = :all,
        xrotation = 45,
        bar_width = 0.5,
        color = 1:size(inv_df, 2),
        xlabel = "Technology",
        ylabel = "Capacity (MW)",
        title = "Generation Investment in VGR",
    )

    Plots.annotate!(
        x,     
        y .+ delta_y,
        cap_str,
        ylim = (0, ymax + 2*delta_y),
        xrotation = 45,
    )

end

a = basic_barchart(vgr_gen_inv)
b = basic_barchart(vgr_sto_inv)

Plots.savefig(a, "gen_1d.png")
Plots.savefig(b, "sto_1d.png")