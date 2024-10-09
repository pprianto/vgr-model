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
1. Model in .mps file

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

    # if model == Model(Gurobi.Optimizer)
    #     compute_conflict!(model)

    #     if get_attribute(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
    #         iis_model, _ = copy_conflict(model)
    #         print(iis_model)
    #     end
    # end
    
    time_solve = time()-start_solve  
    println("---------------------------------------------")
    println("Time needed to solve = $(time_solve / 60) mins")
    println("---------------------------------------------")
    
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
    elpris_df = read_file(joinpath(input_dir, "elpris.csv"))

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
    el_demand_df = read_file(joinpath(input_dir, "el_nodal_demand.csv"))
    heat_demand_df = read_file(joinpath(input_dir, "heat_nodal_demand.csv"))
    h2_demand_df = read_file(joinpath(input_dir, "h2_nodal_demand.csv"))

    # RE profile
    PV_fix_profile = read_file(joinpath(input_dir, "nodal_profile_pv_fixed_$(options.profile_year).csv"))          # fixed axis pv
    PV_opt_profile = read_file(joinpath(input_dir, "nodal_profile_pv_double_axis_$(options.profile_year).csv"))    # opt tracking pv
    WT_on_profile = read_file(joinpath(input_dir, "nodal_profile_onshore_wind_$(options.profile_year).csv"))       # onshore wind
    WT_off_profile = read_file(joinpath(input_dir, "nodal_profile_offshore_wind_$(options.profile_year).csv"))     # offshore wind

    # RE profile in axis arrays
    PVfix = df_to_axisarrays(PV_fix_profile)          # fixed axis pv
    PVopt = df_to_axisarrays(PV_opt_profile)    # opt tracking pv
    WTon = df_to_axisarrays(WT_on_profile)       # onshore wind
    WToff = df_to_axisarrays(WT_off_profile)     # offshore wind

    # collections of input data in NamedTuple
    # price = (; SE3=Array(elpris_df.SE3), NO1=Array(elpris_df.NO1), DK1=Array(elpris_df.DK1))
    # grid_infra = (; subs=substations_df, lines=lines_df, pp=pp_df)
    # tech_props = (; gen=gen_tech_df, sto=sto_tech_df)
    # demand = (; el=el_demand_df, heat=heat_demand_df, h2=h2_demand_df)
    # profiles = (; PVfix=PV_fix_profile, PVopt=PV_opt_profile, WTon=WT_on_profile, WToff=WT_off_profile)

    # collections of input data in Struct
    price = Prices(Array(elpris_df.SE3), Array(elpris_df.NO1), Array(elpris_df.DK1))
    grid_infra = GridInfrastructures(substations_df, lines_df, pp_df)
    tech_props = TechProps(gen_tech_df, sto_tech_df)
    demand = Demands(el_demand_df, heat_demand_df, h2_demand_df)
    # profiles = Profiles(PV_fix_profile, PV_opt_profile, WT_on_profile, WT_off_profile)
    profiles = (; PVfix, PVopt, WTon, WToff)

    return price, grid_infra, tech_props, demand, profiles

end


function set_solver(model, solver::Symbol = :gurobi)

    if solver == :gurobi
        # Workaround for bug:
        # https://discourse.julialang.org/t/how-can-i-clear-solver-attributes-in-jump/57953
        optimizer = optimizer_with_attributes(
                Gurobi.Optimizer,
                "Threads" => 8,                     # shutoff for now, memory limit?
                "BarHomogeneous" => 1,              # 1: enabled
                "Crossover" => 0,                  # 0: disabled
                # "CrossoverBasis" => 0,                  # 0: disabled
                # "BarConvTol" => 1e-6,                  # 0: disabled
                "Method" => 2,                     # -1: auto, 1: dual simplex, 2: barrier
                "Presolve" => 2,                    # 2: aggressive
                "PreSparsify" => 2,                    # 2: aggressive
                "NumericFocus" => 2,
                "NodefileStart" => 0.5,
                "Aggregate" => 1,
                "MemLimit" => 250,
                "ScaleFlag" => 1,
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
                # "Crossover" => 1,       # 0=no, 1=yes
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
    @unpack NODES, 
            TRANSMISSION_NODES,
            COAST_NODES, 
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
            FLEX_TH = sets

    ## Variables

    @unpack existing_generation,
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

    #=------------------------------------------------------------------------------
    EXTRACT OPTIMUM SOLUTION VARIABLES
    ------------------------------------------------------------------------------=#

    if termination_status(model) == MOI.OPTIMAL
        println("---------------------------------------------")
        println("The solver termination status is $(termination_status(model))")
        println("Optimal solution found")
        
        cost = objective_value(model)
        println("The system cost is $(cost) k€.")

        # Investment solutions
        Generation_Investment = DataFrame(NODE = NODES)
        for tech ∈ GEN_TECHS
            Generation_Investment[!, Symbol(tech)] = [value(generation_investment[node, tech]) for node ∈ NODES]
        end

        Storage_Investment = DataFrame(NODE = NODES)
        for tech ∈ STO_TECHS
            Storage_Investment[!, Symbol(tech)] = [value(storage_investment[node, tech]) for node ∈ NODES]
        end

        Generation_Investment = sum_capacities(Generation_Investment)
        Storage_Investment = sum_capacities(Storage_Investment)

        # Dispatch solutions
        Generation_Dispatch = Dict{Symbol, DataFrame}()
        Reactive_Dispatch = Dict{Symbol, DataFrame}()
        Storage_Charge = Dict{Symbol, DataFrame}()
        Storage_Discharge = Dict{Symbol, DataFrame}()
        Storage_Level = Dict{Symbol, DataFrame}()

        for node ∈ NODES
            Generation_Dispatch[node] = DataFrame(
                [value(active_generation[t, node, tech]) for t ∈ PERIODS, tech ∈ GEN_TECHS], 
                GEN_TECHS
            )

            Reactive_Dispatch[node] = DataFrame(
                [value(reactive_generation[t, node, tech]) for t ∈ PERIODS, tech ∈ EL_GEN], 
                EL_GEN
            )

            Storage_Charge[node] = DataFrame(
                [value(storage_charge[t, node, tech]) for t ∈ PERIODS, tech ∈ STO_TECHS], 
                STO_TECHS
            )

            Storage_Discharge[node] = DataFrame(
                [value(storage_discharge[t, node, tech]) for t ∈ PERIODS, tech ∈ STO_TECHS], 
                STO_TECHS
            )

            Storage_Level[node] = DataFrame(
                [value(storage_discharge[t, node, tech]) for t ∈ PERIODS, tech ∈ STO_TECHS], 
                STO_TECHS
            )
        end

        # Long DataFrame format for dispatch
        # for consideration
        # to adjust accordingly

        combined_df = DataFrame(Node=Symbol[], Tech=Symbol[], Period=Int[], Value=Float64[])

        for (node, df) in Generation_Dispatch
            for (tech_idx, tech) in enumerate(GEN_TECHS)
                for period in PERIODS
                    push!(combined_df, (Node=node, Tech=tech, Period=period,Value=df[period, tech_idx]))
                end
            end
        end

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
                    # combined_df
        )
        
        return results

    else
        println("No optimal solution found")
        
        return nothing
    
    end

end

