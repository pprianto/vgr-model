#=------------------------------------------------------------------------------
------------------------------ UTILITIES ---------------------------------------

Anything to support the model,
such as:
1. get the first sheet name to read the XLSX
2. save file of the results
3. etc...

------------------------------------------------------------------------------=#

using DataFrames, AxisArrays, CSV, XLSX, UnPack, JuMP, HiGHS


function get_df(
    data
)

#=------------------------------------------------------------------------------
Specific for xls file type
------------------------------------------------------------------------------=#

    df = XLSX.readxlsx(data)
    sheet_name = XLSX.sheetnames(df)
    first_sheet = sheet_name[1]
    
    # df_data = DataFrame(XLSX.readtable(data, first_sheet))

    # may cause some type problems -> to be solved later
    df_data = DataFrame(XLSX.readtable(data, first_sheet, infer_eltypes=true)) 

    return df_data

end


function read_file(
    file_path::String
)

#=------------------------------------------------------------------------------
Read either xls or csv file as dataframe
------------------------------------------------------------------------------=#


    if endswith(file_path, ".xlsx") || endswith(file_path, ".xls")
        # Read XLSX file
        df = get_df(file_path)
    
    elseif endswith(file_path, ".csv")
        # Read CSV file
        df = CSV.read(file_path, DataFrame)
    
    else
        error("Unsupported file format. Please provide XLS or CSV file.")
    
    end
    
    return df

end


function carrier_subsets(
    tech_df::DataFrame
)

#=------------------------------------------------------------------------------
filter df based on tech carrier
------------------------------------------------------------------------------=#

    if any(col -> any(x -> contains(string(x), "EL"), col), eachcol(tech_df))
        el_tech = filter(row -> occursin("EL", row.Carrier), tech_df)
        el_subset = el_tech.Tech
    end

    if any(col -> any(x -> contains(string(x), "HEAT"), col), eachcol(tech_df))
        heat_tech = filter(row -> occursin("HEAT", row.Carrier), tech_df)
        heat_subset = heat_tech.Tech
    end

    if any(col -> any(x -> contains(string(x), "H2"), col), eachcol(tech_df))    
        h2_tech = filter(row -> occursin("H2", row.Carrier), tech_df)
        H2_subset = h2_tech.Tech
    end

    return el_subset, heat_subset, H2_subset

end


function rename_pp(
    pp_df::DataFrame
)
    
    tech_names = Dict(
        "wind" => "WON",            # change into onshore wind
        "hydro" => "HYD",           # hydro
        "waste;biomass" => "WCHP",  # waste CHP
        "biomass" => "WCCHP",       # wood chips biomass chp
        "electric" => "EB",         # electric boiler
        "biofuel" => "GEBG",        # gas engines biogas
        "solar" => "PVUTIL",        # utility scale PV
        "gas;oil" => "CCGT",        # combined cycle gas turbine
        "gas;wood;oil" => "CCGT",   # combined cycle gas turbine
        "oil" => "GTSC"             # gas turbine simple cycle
    )

    pp_df[!, :tech] = [tech_names[name] for name in pp_df[:, :tech]]


    # replace missing values in installed capacity column with 0.0
    power_cols = [:el_MW, :heat_MW]

    for col in power_cols
        pp_df[!, col] = replace(pp_df[:, col], missing => 0.0)
    end

    return pp_df

end


function admittance_matrix(
    lines_df::DataFrame
)

#=------------------------------------------------------------------------------
Function to create admittance matrix of the network
adapted from the JuMP documentation
https://jump.dev/JuMP.jl/stable/tutorials/applications/optimal_power_flow/
------------------------------------------------------------------------------=#

    # get number of nodes and lines
    max_node_num = max(
                        maximum(lines_df.station_from), 
                        maximum(lines_df.station_to)
                    )

    min_node_num = min(
                        minimum(lines_df.station_from), 
                        minimum(lines_df.station_to)
                    )

    no_nodes = max_node_num
    no_lines = size(lines_df, 1)
    
    # arrays for row/column
    nodes = [i for i in min_node_num : max_node_num]


    # construct incidence matrix
    incidence_matrix = SparseArrays.sparse(
                            lines_df.station_from,      # vertice
                            1:no_lines,                 # edge
                            1,                          # to
                            no_nodes,                   # row size
                            no_lines                    # from size
                            ) +
                        SparseArrays.sparse(
                            lines_df.station_to,        # vertice
                            1:no_lines,                 # edge
                            -1,                         # from
                            no_nodes,                   # row size
                            no_lines                    # from size
                        )

    # line impedance data in complex
    z_total = lines_df.r_total .+ im * lines_df.x_total

    # admittance matrix, change from SparseArrays to Array
    Ybus = incidence_matrix * SparseArrays.spdiagm(1 ./ z_total) * incidence_matrix'
    Ybus = Array(Ybus)

    # conductance and susceptance matrix    
    G = real(Ybus)
    B = imag(Ybus)

    # change to dataframe for later uses (parameters in the model)
    Ybus = DataFrame(Ybus, Symbol.(nodes))
    G = DataFrame(G, Symbol.(nodes))
    B = DataFrame(B, Symbol.(nodes))

    return Ybus, G, B

end