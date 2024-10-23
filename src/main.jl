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


Example of run:

@time m = run_model(:copt);

@time results = query_solutions(
    m.model,
    m.sets,
    m.vars
);

vgr_gen_cap = vgr_investment(results.Generation_Capacity)
vgr_gen_inv = vgr_investment(results.Generation_Investment)
vgr_sto_inv = vgr_investment(results.Storage_Investment)

a = basic_barchart(vgr_gen_cap)
b = basic_barchart(vgr_gen_inv)
c = basic_barchart(vgr_sto_inv)

Plots.savefig(a, "gen_1d.png")
Plots.savefig(b, "sto_1d.png")

------------------------------------------------------------------------------=#

using DataFrames, CSV, XLSX, Printf, JLD2
using JuMP, HiGHS, Ipopt, Gurobi, COPT
import AxisArrays, SparseArrays, LinearAlgebra
import Plots, CairoMakie
CairoMakie.activate!()

include("structs.jl")
include("utilities.jl")
include("params_vars.jl")
include("constraints.jl")
include("model.jl")

# Set directories
const current_dir::String = pwd()
mkpath("modelinput")                                            # input folder, should already be created during retrieve OSM process
mkpath("results")                                               # results folder
const input_dir::String = joinpath(current_dir, "modelinput")
const results_dir::String = joinpath(current_dir, "results")

# Model options
# decided as global for now so that can be called in functions
const options::ModelOptions = ModelOptions(
                                            run = :trial,           # :full or :trial :test
                                            CO2_budget = 0.0,       # in kgCO2
                                            target_year = 2050,     # model target year
                                            profile_year = 2019,    # model profile
                                            el_price_year = 1991,   # future el price model year (1991 or 1992)
                                            scenario = :beta,       # :now, :alpha, :beta, with beta as highest demand increase
                                            power_flow = :lac,      # linearised AC (LAC) or DC
                                            Discount_rate = 0.05,   # assumed discount rate
                                            El_Heat_Tax = 60.0,     # assumed tax to generate electricity, on top of el price (€/MWh)
                                            CO2_limit = :fee,       # constraint to limit CO2 emitting technology (:ton or :fee)
                                            Emission_fee = 150.0,   # assumed CO2 emission fee (€/tonCO2)
                                            Vnom = 130.0,           # system voltage level (kV)
                                            FlexLim = :yes,         # option if flexlim constraint is included
                                            EV = :yes,              # option if EV is included
)          

if options.EV == :yes
    const EV_options::EVOptions = EVOptions()
end

nothing


