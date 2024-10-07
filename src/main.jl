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
using JuMP, HiGHS, Ipopt, Gurobi, COPT
import AxisArrays, SparseArrays, LinearAlgebra
import Plots, CairoMakie
CairoMakie.activate!()

include("structs.jl")
include("utilities.jl")
include("params_vars.jl")
include("constraints.jl")
include("model.jl")

# Set directory

const current_dir::String = pwd()
mkpath("modelinput")                                            # input folder, should already be created during retrieve OSM process
mkpath("results")                                               # results folder
const input_dir::String = joinpath(current_dir, "modelinput")
const results_dir::String = joinpath(current_dir, "results")
const options::ModelOptions = ModelOptions(run=:trial, FlexLim=:yes, EV=:yes)          # decided as global for now so that can be called in functions
if options.EV == :yes
    const EV::EVOptions = EVOptions()          # decided as global for now so that can be called in functions
end

@time m = run_model(:copt);


(; model, sets, vars) = m

@time results = query_solutions(
    model,
    sets,
    vars
);


vgr_gen_inv = vgr_investment(results.Generation_Investment)
vgr_sto_inv = vgr_investment(results.Storage_Investment)

a = basic_barchart(vgr_gen_inv)
b = basic_barchart(vgr_sto_inv)

# Plots.savefig(a, "gen_1d.png")
# Plots.savefig(b, "sto_1d.png")
