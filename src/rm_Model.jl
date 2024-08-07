function run_Model(
    model::Model,
    price::NamedTuple,
    grid_infra::NamedTuple,
    tech_props::NamedTuple,
    demand::NamedTuple,
    profiles::NamedTuple   
)
#=------------------------------------------------------------------------------
-------------------------------- MODEL -----------------------------------------

Define core Model
Taking the following:
1. Sets
2. Params
3. Vars
4. Constraints
5. Query results

Return:
1. Model
2. Results in NamedTuple

------------------------------------------------------------------------------=#
    println("---------------------------------------------")
    println("Read input data")
    println("---------------------------------------------")
    start_sets = time()

    sets, params = make_Params_Sets(
        price,
        grid_infra,
        tech_props,
        demand
    )

    time_sets = time()-start_sets  

    println("---------------------------------------------")
    println("Time needed to read input = $(time_sets / 60) mins")
    println("---------------------------------------------")

    #------------------------------------------------------------------------------=#
    println("---------------------------------------------")
    println("Define Variables")
    println("---------------------------------------------")
    start_vars = time()

    vars = make_Variables(
        model,
        sets,
        params
    )

    time_vars = time()-start_vars  
    println("---------------------------------------------")
    println("Time needed to define variables = $(time_vars / 60) mins")
    println("---------------------------------------------")

    #------------------------------------------------------------------------------=#
    println("---------------------------------------------")
    println("Define Constraints")
    println("---------------------------------------------")
    start_const = time()

    constraints = make_Constraints(
        model,
        sets,
        params,
        vars,
        grid_infra,
        profiles
    )

    time_consts = time()-start_const  
    println("---------------------------------------------")
    println("Time needed to define constraints = $(time_consts / 60) mins")
    println("---------------------------------------------")

    #------------------------------------------------------------------------------=#    
    @unpack total_cost = vars

    @objective model Min begin
        total_cost
    end

    println("---------------------------------------------")
    println("Start Solving")
    println("---------------------------------------------")
    start_solve = time()

    optimize!(model)

    if model == Model(Gurobi.Optimizer)
        compute_conflict!(model)

        if get_attribute(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
            iis_model, _ = copy_conflict(model)
            print(iis_model)
        end
    end

    println("---------------------------------------------")
    println("The solver termination status is $(termination_status(model))")
    
    cost = objective_value(model)
    println("The cost of system is $(cost) M€.")
    
    time_solve = time()-start_solve  
    println("---------------------------------------------")
    println("Time needed to solve = $(time_solve / 60) mins")
    println("---------------------------------------------")
    
    #=------------------------------------------------------------------------------
    RETURN
    ------------------------------------------------------------------------------=#  

    results = query_solutions(
        model,
        sets,
        vars
    )

    times = (;
            time_sets,
            time_vars,
            time_consts,
            time_solve
    )

    return model, results, times

end


function query_solutions(
    model::Model,
    sets::NamedTuple,
    vars::NamedTuple
)
#=------------------------------------------------------------------------------
---------------------------- QUERY SOLUTIONS -----------------------------------

Extract the optimum solutions
Taking the following:
1. Model structure
2. Model Sets
3. Model Variables

Return:
optimal solution in Dict or DataFrame, TODO: to consider other data structure?
which then converted into NamedTuple
including:
1. Investments
2. Dispatch
3. Export/Import
4. Nodal voltage
5. Power flow in lines

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
            ARCS_FR,
            ARCS_TO,
            CHP,
            FC,
            WIND,
            PV,
            HP,
            BOILER,
            EC,
            FLEX_TH = sets

    ## Variables

    @unpack total_cost, 
            generation_capacity, 
            storage_capacity, 
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
        println("Optimal solution found")

        # Investment solutions
        Generation_Capacity = Dict{String, DataFrame}()
        Storage_Capacity = Dict{String, DataFrame}()

        for (_, node) in enumerate(NODES)
            Generation_Capacity[node] = DataFrame(
                [Symbol(tech) => value(generation_capacity[node, tech]) for tech in GEN_TECHS]
            )

            Storage_Capacity[node] = DataFrame(
                [Symbol(tech) => value(storage_capacity[node, tech]) for tech in STO_TECHS]
            )
        end

        # Dispatch solutions
        Generation_Dispatch = Dict{String, DataFrame}()
        Storage_Charge = Dict{String, DataFrame}()
        Storage_Discharge = Dict{String, DataFrame}()
        Storage_Level = Dict{String, DataFrame}()

        for (_, node) in enumerate(NODES)
            Generation_Dispatch[node] = DataFrame(
                [value(active_generation[node, tech, t]) for t in PERIODS, tech in GEN_TECHS], 
                GEN_TECHS
            )

            Storage_Charge[node] = DataFrame(
                [value(storage_charge[node, tech, t]) for t in PERIODS, tech in STO_TECHS], 
                STO_TECHS
            )

            Storage_Discharge[node] = DataFrame(
                [value(storage_discharge[node, tech, t]) for t in PERIODS, tech in STO_TECHS], 
                STO_TECHS
            )

            Storage_Level[node] = DataFrame(
                [value(storage_discharge[node, tech, t]) for t in PERIODS, tech in STO_TECHS], 
                STO_TECHS
            )
        end

        # Import Export solutions

        Export_to = DataFrame()
        Import_from = DataFrame()
        Export_Import = DataFrame()

        for (_, node) in enumerate(TRANSMISSION_NODES)
            Export_to[!, Symbol(node)] = [value(export_to[node, t]) for t in PERIODS]
            Import_from[!, Symbol(node)] = [value(import_from[node, t]) for t in PERIODS]
            Export_Import[!, Symbol(node)] = [value(import_export[node, t]) for t in PERIODS]
        end

        # Voltage
        Nodal_Voltage = DataFrame()

        for (_, node) in enumerate(NODES)
            Nodal_Voltage[!, Symbol(node)] = [value(nodal_voltage[node, t]) for t in PERIODS]
        end

        # Power Flow
        Active_Flow = DataFrame()

        for (_, line) in enumerate(LINES)
            Active_Flow[!, Symbol(line)] = [value(active_flow[line, t]) for t in PERIODS]
        end

        results =   (; 
                    Generation_Capacity, 
                    Storage_Capacity,
                    Generation_Dispatch,
                    Storage_Charge,
                    Storage_Discharge,
                    Storage_Level,
                    Export_to,
                    Import_from,
                    Export_Import,
                    Nodal_Voltage,
                    Active_Flow
        )
        
        
        return results

    else
        println("No optimal solution found")
        
        return nothing
    
    end

end