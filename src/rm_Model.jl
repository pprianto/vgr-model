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

    vars = make_Constraints(
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
    
    time_solve = time()-start_solve  
    println("---------------------------------------------")
    println("Time needed to solve = $(time_solve / 60) mins")
    println("---------------------------------------------")
    
    #=------------------------------------------------------------------------------
    RETURN
    ------------------------------------------------------------------------------=#  

    # results = query_solutions(
    #     model,
    #     sets,
    #     vars
    # )

    times = (;
            time_sets,
            time_vars,
            time_consts,
            time_solve
    )

    model_struct = Model_Struct(
                                model,
                                sets,
                                params,
                                vars,
                                times
    )
    
    return model_struct

    # return model, results, times

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
        Generation_Dispatch = Dict{String, DataFrame}()
        Reactive_Dispatch = Dict{String, DataFrame}()
        Storage_Charge = Dict{String, DataFrame}()
        Storage_Discharge = Dict{String, DataFrame}()
        Storage_Level = Dict{String, DataFrame}()

        for node ∈ NODES
            Generation_Dispatch[node] = DataFrame(
                [value(active_generation[node, tech, t]) for t ∈ PERIODS, tech ∈ GEN_TECHS], 
                GEN_TECHS
            )

            Reactive_Dispatch[node] = DataFrame(
                [value(reactive_generation[node, tech, t]) for t ∈ PERIODS, tech ∈ EL_GEN], 
                EL_GEN
            )

            Storage_Charge[node] = DataFrame(
                [value(storage_charge[node, tech, t]) for t ∈ PERIODS, tech ∈ STO_TECHS], 
                STO_TECHS
            )

            Storage_Discharge[node] = DataFrame(
                [value(storage_discharge[node, tech, t]) for t ∈ PERIODS, tech ∈ STO_TECHS], 
                STO_TECHS
            )

            Storage_Level[node] = DataFrame(
                [value(storage_discharge[node, tech, t]) for t ∈ PERIODS, tech ∈ STO_TECHS], 
                STO_TECHS
            )
        end

        # Long DataFrame format for dispatch
        # for consideration
        # to adjust accordingly

        combined_df = DataFrame(Node=String[], Tech=String[], Period=Int[], Value=Float64[])

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
            Export_to[!, Symbol(node)] = [value(export_to[node, t]) for t ∈ PERIODS]
            Import_from[!, Symbol(node)] = [value(import_from[node, t]) for t ∈ PERIODS]
            Export_Import[!, Symbol(node)] = [value(import_export[node, t]) for t ∈ PERIODS]
        end

        # Voltage
        Nodal_Voltage = DataFrame()
        Nodal_Angle = DataFrame()

        for node ∈ NODES
            Nodal_Voltage[!, Symbol(node)] = [value(nodal_voltage[node, t]) for t ∈ PERIODS]
            Nodal_Angle[!, Symbol(node)] = [value(nodal_angle[node, t]) for t ∈ PERIODS]
        end

        # Power Flow
        Active_Flow = DataFrame()
        Reactive_Flow = DataFrame()

        for line ∈ LINES
            Active_Flow[!, Symbol(line)] = [value(active_flow[line, t]) for t ∈ PERIODS]
            Reactive_Flow[!, Symbol(line)] = [value(reactive_flow[line, t]) for t ∈ PERIODS]
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

