function run_model(
    # model::Model,
    # price::NamedTuple,
    # grid_infra::NamedTuple,
    # tech_props::NamedTuple,
    # demand::NamedTuple,
    # profiles::NamedTuple,
    # options=ModelOptions(),
    solver::Symbol = :gurobi,
)
#=------------------------------------------------------------------------------
-------------------------------- MODEL -----------------------------------------

Define the problem Model
Taking the following:
1. Sets
2. Params
3. Vars
4. Constraints

Return:
1. Model in NamedTuple

------------------------------------------------------------------------------=#

    model = Model()

    #=---------------------------------------------
    INPUT DATA
    ---------------------------------------------=#
    println("---------------------------------------------")
    println("Read input data")
    println("---------------------------------------------")

    price, grid_infra, tech_props, demand, profiles = read_input_data()

    start_sets = time()

    @time sets, params = make_sets(price, grid_infra, tech_props, demand)

    time_sets = time()-start_sets  

    println("---------------------------------------------")
    println("Time needed to read input = $(time_sets / 60) mins")
    println("---------------------------------------------")

    #------------------------------------------------------------------------------=#
    println("---------------------------------------------")
    println("Define Variables")
    println("---------------------------------------------")
    start_vars = time()

    @time vars = make_variables(model, sets, params)

    time_vars = time()-start_vars  
    println("---------------------------------------------")
    println("Time needed to define variables = $(time_vars / 60) mins")
    println("---------------------------------------------")

    #------------------------------------------------------------------------------=#
    println("---------------------------------------------")
    println("Define Constraints")
    println("---------------------------------------------")
    start_const = time()

    @time gen_constraints(model, sets, params, vars, grid_infra, profiles)
    if options.FlexLim == :yes
        @time flex_lim_constraints(model, sets, params, vars, grid_infra, profiles)
    end
    @time sto_constraints(model, sets, params, vars, grid_infra, profiles)
    @time enbal_constraints(model, sets, params, vars, grid_infra, profiles)
    @time power_flow_constraints(model, sets, params, vars, grid_infra, profiles)
    @time cost_constraints(model, sets, params, vars, grid_infra, profiles)

    time_consts = time()-start_const  
    println("---------------------------------------------")
    println("Time needed to define constraints = $(time_consts / 60) mins")
    println("---------------------------------------------")

    # ------------------------------------------------------------------------------=#    

    println("---------------------------------------------")
    println("Start Solving")
    println("---------------------------------------------")
    start_solve = time()

    set_solver(model, solver)

    optimize!(model)

    time_solve = time()-start_solve  

    if time_solve ≥ 3600
        println("---------------------------------------------")
        println("Time needed to solve = $(time_solve / 3600) hours")
        println("---------------------------------------------")
    else
        println("---------------------------------------------")
        println("Time needed to solve = $(time_solve / 60) mins")
        println("---------------------------------------------")
    end    
    #=------------------------------------------------------------------------------
    RETURN
    ------------------------------------------------------------------------------=#  

    times = (;
            time_sets,
            time_vars,
            time_consts,
            time_solve
    )

    # model_struct = ModelStruct(model, sets, params, vars, times)
    model_struct = (; model, sets, params, vars, times)
    
    return model_struct

end


function read_input_data(
    # options=ModelOptions()
)
    # electricity price
    elpris_df = read_file(joinpath(input_dir, "elpris_$(options.el_price_year).csv"))

    # electrical infrastructure
    substations_df = read_file(joinpath(input_dir, "subs_final.csv"))
    lines_df = read_file(joinpath(input_dir, "lines_final.csv"))
    pp_df = read_file(joinpath(input_dir, "pp_final.csv"))
    pp_df = process_pp(pp_df)                                                                    # rename the tech according to the index sets

    # technology properties
    # assume 2050
    gen_tech_df = read_file(joinpath(input_dir, "gen_tech_$(options.target_year).csv"))
    sto_tech_df = read_file(joinpath(input_dir, "sto_tech_$(options.target_year).csv"))

    # demand data
    # currently in hourly period
    el_demand_df = read_file(joinpath(input_dir, "el_nodal_demand_$(options.scenario).csv"))
    el_demand_df = new_industries_el_demand(el_demand_df)
    heat_demand_df = read_file(joinpath(input_dir, "heat_nodal_demand_$(options.scenario).csv"))
    h2_demand_df = read_file(joinpath(input_dir, "h2_nodal_demand_$(options.scenario).csv"))

    # RE profiles
    PV_fix_profile = read_file(joinpath(input_dir, "nodal_profile_pv_fixed_$(options.profile_year).csv"))          # fixed axis pv
    PV_opt_profile = read_file(joinpath(input_dir, "nodal_profile_pv_double_axis_$(options.profile_year).csv"))    # opt tracking pv
    WT_on_profile = read_file(joinpath(input_dir, "nodal_profile_onshore_wind_$(options.profile_year).csv"))       # onshore wind
    WT_off_profile = read_file(joinpath(input_dir, "nodal_profile_offshore_wind_$(options.profile_year).csv"))     # offshore wind

    # RE profile in axis arrays
    PVfix = df_to_axisarrays(PV_fix_profile)    # fixed axis pv
    PVopt = df_to_axisarrays(PV_opt_profile)    # opt tracking pv
    WTon = df_to_axisarrays(WT_on_profile)      # onshore wind
    WToff = df_to_axisarrays(WT_off_profile)    # offshore wind

    # collections of input data for return
    price = Prices(Array(elpris_df.SE3), Array(elpris_df.NO1), Array(elpris_df.DK1))
    grid_infra = GridInfrastructures(substations_df, lines_df, pp_df)
    tech_props = TechProps(gen_tech_df, sto_tech_df)
    demand = Demands(el_demand_df, heat_demand_df, h2_demand_df)
    profiles = (; PVfix, PVopt, WTon, WToff)
    
    return price, grid_infra, tech_props, demand, profiles

end


function set_solver(model, solver::Symbol = :gurobi)

    if solver == :gurobi
        # Workaround for bug:
        # https://discourse.julialang.org/t/how-can-i-clear-solver-attributes-in-jump/57953
        optimizer = optimizer_with_attributes(
                Gurobi.Optimizer,
                "Threads" => 24,                     # shutoff for now, memory limit?
                "BarHomogeneous" => 1,              # 1: enabled
                "Crossover" => 0,                  # 0: disabled
                # "CrossoverBasis" => 0,                  # 0: disabled
                # "BarConvTol" => 1e-6,                  # 0: disabled
                "Method" => 2,                     # -1: auto, 1: dual simplex, 2: barrier
                "Presolve" => 2,                    # 2: aggressive
                "PreSparsify" => 2,                    # 2: aggressive
                "NumericFocus" => 1,
                # "NodefileStart" => 0.5,
                "Aggregate" => 2,
                # "MemLimit" => 250,
                "ScaleFlag" => 3,
                "DisplayInterval" => 300,
        )
        
        # set_silent(model)
        # log_file = joinpath(results_dir, "1year_log.txt")
        # set_optimizer_attribute(model, "LogFile", log_file)

    elseif solver == :copt
        optimizer = optimizer_with_attributes(
                COPT.Optimizer,
                "LpMethod" => 2,        # -1=simple auto, 1=Dual simplex, 2=Barrier, 3=Crossover, 4=Concurrent, 5=heuristic auto, 6=PDLP
                "BarIterLimit" => 1e9,
                # "BarHomogeneous" => 1,  # 0=no, 1=yes
                # "BarOrder" => 1,       # -1=auto, 0=Approximate Minimum Degree, 1=Nested Dissection
                # "BarStart" => 2,       # -1=auto, 0=Asimple, 1=Mehrotra, 2=Modified Mehrotra
                "Crossover" => 0,       # 0=no, 1=yes
                # "Presolve" => 4,        # -1=auto, 0=off, 1=fast, 2=normal, 3=aggressive, 4=unlimited (until nothing else possible)
                # "Scaling" => 1,         # -1=auto, 0=no, 1=yes
                # "GPUMode" => 1,       # 0=CPU, 1=GPU (used with PDLP algoritm, only for machine 41)
                # "PDLPTol" => 1e-7,
                "Threads" => 24,
                "BarThreads" => 24,
                "SimplexThreads" => 24,
        )

    elseif solver == :highs
        optimizer = optimizer_with_attributes(
            HiGHS.Optimizer,
            # "threads" => 0,
            "presolve" => "on",
            "solver" => "ipm",          # simplex, ipm=barrier, pdlp
            "run_crossover" => "off",
        )

    else
        @error "No solver named $(solver)."
    end

    set_optimizer(model, optimizer)

end


function query_solutions(
    model::Model,
    # sets::ModelSets,
    # vars::ModelVariables,
    sets::NamedTuple,
    vars::NamedTuple,    
)
#=------------------------------------------------------------------------------
---------------------------- QUERY SOLUTIONS -----------------------------------

Extract the optimum solutions
Taking the following:
1. Model structure
2. Model sets
3. Model variables

Return:
optimal solution in Dict or DataFrame, TODO: to consider other data structure?
which then converted into NamedTuple
including:
1. Investments
2. Dispatch (P&Q)
3. Export/Import
4. Voltage magnitude and angle
5. Power flow in lines (P&Q)

Considerations for data structure:
1. easier to save as csv or xls
2. easier for plotting

------------------------------------------------------------------------------=#

    #=------------------------------------------------------------------------------
    MODEL SETS & VARIABLES
    ------------------------------------------------------------------------------=#
    (;  NODES, 
        TRANSMISSION_NODES,
        GEN_TECHS, 
        EL_GEN, 
        STO_TECHS, 
        PERIODS, 
        LINES,
        STO_EN 
    ) = sets

    ## Variables

    (;  total_cost,
        capex,
        fix_om,
        fuel_cost,
        var_om,
        start_part_costs,
        exp_imp_costs,
        tax_cost,
        CO2_cost,
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
        reactive_flow
    ) = vars

    #=------------------------------------------------------------------------------
    EXTRACT OPTIMUM SOLUTION VARIABLES
    ------------------------------------------------------------------------------=#

    if termination_status(model) == MOI.OPTIMAL
        println("---------------------------------------------")
        println("The solver termination status is $(termination_status(model))")
        println("Optimal solution found")
        
        cost = objective_value(model)
        println("The system cost is $(cost) k€.")

        # Cost variables
        Costs = DataFrame(
            :System_cost => value(total_cost),
            :CAPEX_cost => value(capex),
            :Fix_OM_cost => value(fix_om),
            :Fuel_cost => value(fuel_cost),
            :Var_OM_cost => value(var_om),
            :Start_Part_cost => value(start_part_costs),
            :Export_Import_cost => value(exp_imp_costs),
            :Taxes_cost => value(tax_cost),
            :CO2_cost => value(CO2_cost)
        )

        # Investment solutions
        Generation_Capacity = DataFrame(NODE = NODES)
        for tech ∈ GEN_TECHS
            Generation_Capacity[!, Symbol(tech)] = [value((existing_generation[node, tech] + generation_investment[node, tech])) for node ∈ NODES]
        end

        Generation_Investment = DataFrame(NODE = NODES)
        for tech ∈ GEN_TECHS
            Generation_Investment[!, Symbol(tech)] = [value(generation_investment[node, tech]) for node ∈ NODES]
        end

        Storage_Investment = DataFrame(NODE = NODES)
        for tech ∈ STO_TECHS
            Storage_Investment[!, Symbol(tech)] = [value(storage_investment[node, tech]) for node ∈ NODES]
        end

        Generation_Capacity = sum_capacities(Generation_Capacity)
        Generation_Investment = sum_capacities(Generation_Investment)
        Storage_Investment = sum_capacities(Storage_Investment)

        # Dispatch solutions
        Generation_Dispatch = DataFrame(Node=Symbol[], Tech=Symbol[], Period=Int[], Dispatch=Float64[])
        Reactive_Dispatch = DataFrame(Node=Symbol[], Tech=Symbol[], Period=Int[], Dispatch=Float64[])
        Storage_Charge = DataFrame(Node=Symbol[], Tech=Symbol[], Period=Int[], Charge=Float64[])
        Storage_Discharge = DataFrame(Node=Symbol[], Tech=Symbol[], Period=Int[], Discharge=Float64[])
        Storage_Level = DataFrame(Node=Symbol[], Tech=Symbol[], Period=Int[], Level=Float64[])


        for node ∈ NODES
            for tech ∈ GEN_TECHS  
                for time ∈ PERIODS
                    push!(Generation_Dispatch, (node, tech, time, value(active_generation[time, node, tech])))
                end
            end
        end

        # aggregate the nodes dispatch into total (VGR)
        VGR_gen_dispatch = combine(groupby(Generation_Dispatch, [:Tech, :Period]), :Dispatch => sum => :Dispatch)
        VGR_gen_dispatch[!, :Node] .= :VGR
        VGR_gen_dispatch = select(VGR_gen_dispatch, :Node, :Tech, :Period, :Dispatch)
        Generation_Dispatch = vcat(Generation_Dispatch, VGR_gen_dispatch)

        for node ∈ NODES
            for tech ∈ EL_GEN  
                for time ∈ PERIODS
                    push!(Reactive_Dispatch, (node, tech, time, value(reactive_generation[time, node, tech])))
                end
            end
        end

        # aggregate the nodes dispatch into total (VGR)
        VGR_rea_dispatch = combine(groupby(Reactive_Dispatch, [:Tech, :Period]), :Dispatch => sum => :Dispatch)
        VGR_rea_dispatch[!, :Node] .= :VGR
        VGR_rea_dispatch = select(VGR_rea_dispatch, :Node, :Tech, :Period, :Dispatch)
        Reactive_Dispatch = vcat(Reactive_Dispatch, VGR_rea_dispatch)

        for node ∈ NODES
            for tech ∈ STO_EN  
                for time ∈ PERIODS
                    push!(Storage_Charge, (node, tech, time, value(storage_charge[time, node, tech])))
                    push!(Storage_Discharge, (node, tech, time, value(storage_discharge[time, node, tech])))
                    push!(Storage_Level, (node, tech, time, value(storage_level[time, node, tech])))

                end
            end
        end

        # aggregate the nodes charge, discharge, level into total (VGR)
        VGR_sto_charge = combine(groupby(Storage_Charge, [:Tech, :Period]), :Charge => sum => :Charge)
        VGR_sto_charge[!, :Node] .= :VGR
        VGR_sto_charge = select(VGR_sto_charge, :Node, :Tech, :Period, :Charge)
        Storage_Charge = vcat(Storage_Charge, VGR_sto_charge)

        VGR_sto_discharge = combine(groupby(Storage_Discharge, [:Tech, :Period]), :Discharge => sum => :Discharge)
        VGR_sto_discharge[!, :Node] .= :VGR
        VGR_sto_discharge = select(VGR_sto_discharge, :Node, :Tech, :Period, :Discharge)
        Storage_Discharge = vcat(Storage_Discharge, VGR_sto_discharge)

        VGR_sto_lv = combine(groupby(Storage_Level, [:Tech, :Period]), :Level => sum => :Level)
        VGR_sto_lv[!, :Node] .= :VGR
        VGR_sto_lv = select(VGR_sto_lv, :Node, :Tech, :Period, :Level)
        Storage_Level = vcat(Storage_Level, VGR_sto_lv)

        # Import Export solutions
        Export_to = DataFrame()
        Import_from = DataFrame()
        Export_Import = DataFrame()

        for node ∈ TRANSMISSION_NODES
            Export_to[!, Symbol(node)] = [value(export_to[t, node]) for t ∈ PERIODS]
            Import_from[!, Symbol(node)] = [value(import_from[t, node]) for t ∈ PERIODS]
            Export_Import[!, Symbol(node)] = [value(import_export[t, node]) for t ∈ PERIODS]
        end

        # Voltage
        Nodal_Voltage = DataFrame()
        Nodal_Angle = DataFrame()

        for node ∈ NODES
            Nodal_Voltage[!, Symbol(node)] = [value(nodal_voltage[t, node]) for t ∈ PERIODS]
            Nodal_Angle[!, Symbol(node)] = [value(nodal_angle[t, node]) for t ∈ PERIODS]
        end

        # Power Flow
        Active_Flow = DataFrame()
        Reactive_Flow = DataFrame()

        for line ∈ LINES
            Active_Flow[!, Symbol(line)] = [value(active_flow[t, line]) for t ∈ PERIODS]
            Reactive_Flow[!, Symbol(line)] = [value(reactive_flow[t, line]) for t ∈ PERIODS]
        end

        results =   (; 
                    Costs,
                    Generation_Capacity,
                    Generation_Investment, 
                    Storage_Investment,
                    Generation_Dispatch,
                    Reactive_Dispatch,
                    Storage_Charge,
                    Storage_Discharge,
                    Storage_Level,
                    Export_to,
                    Import_from,
                    Export_Import,
                    Nodal_Voltage,
                    Nodal_Angle,
                    Active_Flow,
                    Reactive_Flow
        )
        
        return results

    else
        println("No optimal solution found")
        
        return nothing
    
    end

end


function save_variables(
    results::NamedTuple
)
#=------------------------------------------------------------------------------
----------------------------- SAVE VARIABLES -----------------------------------

Extract the optimum solutions
and save into files
so that can be post-processed further or in another computer

Considerations for data structure:
1. 2-dimensional variables in CSV
2. 3-dimensional variables in JLD2

------------------------------------------------------------------------------=#

    (;  Costs,
        Generation_Capacity,
        Generation_Investment, 
        Storage_Investment,
        Generation_Dispatch,
        Reactive_Dispatch,
        Storage_Charge,
        Storage_Discharge,
        Storage_Level,
        Export_to,
        Import_from,
        Export_Import,
        Nodal_Voltage,
        Nodal_Angle,
        Active_Flow,
        Reactive_Flow
    ) = results

    # Save decision variables
    CSV.write(joinpath(results_dir, "Costs.csv"), Costs; delim=";")                        # 1
    CSV.write(joinpath(results_dir, "Gen_Cap.csv"), Generation_Capacity; delim=";")
    CSV.write(joinpath(results_dir, "Gen_Inv.csv"), Generation_Investment; delim=";")
    CSV.write(joinpath(results_dir, "Sto_Inv.csv"), Storage_Investment; delim=";")
    CSV.write(joinpath(results_dir, "Export.csv"), Export_to; delim=";")
    CSV.write(joinpath(results_dir, "Import.csv"), Import_from; delim=";")
    CSV.write(joinpath(results_dir, "Exp_Imp.csv"), Export_Import; delim=";")
    CSV.write(joinpath(results_dir, "VM_node.csv"), Nodal_Voltage; delim=";")
    CSV.write(joinpath(results_dir, "VA_node.csv"), Nodal_Angle; delim=";")
    CSV.write(joinpath(results_dir, "P_line.csv"), Active_Flow; delim=";")
    CSV.write(joinpath(results_dir, "Q_line.csv"), Reactive_Flow; delim=";")

    CSV.write(joinpath(results_dir, "Gen_Disp.csv"), Generation_Dispatch; delim=";")
    CSV.write(joinpath(results_dir, "Q_Disp.csv"), Reactive_Dispatch; delim=";")
    CSV.write(joinpath(results_dir, "Sto_Ch.csv"), Storage_Charge; delim=";")
    CSV.write(joinpath(results_dir, "Sto_Dch.csv"), Storage_Discharge; delim=";")
    CSV.write(joinpath(results_dir, "Sto_Lv.csv"), Storage_Level; delim=";")

    # jldsave(joinpath(results_dir, "Gen_Disp.jld2"); Generation_Dispatch)
    # jldsave(joinpath(results_dir, "Q_Disp.jld2"); Reactive_Dispatch)
    # jldsave(joinpath(results_dir, "Sto_Ch.jld2"); Storage_Charge)
    # jldsave(joinpath(results_dir, "Sto_Dch.jld2"); Storage_Discharge)
    # jldsave(joinpath(results_dir, "Sto_Lv.jld2"); Storage_Level)                # 16

end