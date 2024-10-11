#=------------------------------------------------------------------------------
------------------------------ UTILITIES ---------------------------------------

Anything to support the model,
such as:
1. get the first sheet name to read the XLSX
2. save file of the results
3. etc...

------------------------------------------------------------------------------=#
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

    string_columns = names(df, col -> occursin("String", string(eltype(df[!, col]))))

    for col in string_columns
        df[!, col] = [ismissing(val) ? missing : String(val) for val in df[!, col]]
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


function process_pp(
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

    # Stenungsund decommissionioned, power set to 0
    stenungsund_pp = findall(x -> x in [15009765], pp_df.pp_id)
    pp_df[stenungsund_pp, :el_MW] .= 0.0

    # hydro capacities adjusted based on
    # https://vattenkraft.info/
    # crosscheck with vattenfall or relevant websites
    hydro_change = findall(x -> x in [
                                    88144184,   # olidans
                                    10884545,   # vargons
                                    528968206,  # kungfors
                                    125589869   # lilla edet
                                    ],
                            pp_df.pp_id
    )

    # ordered based on appearance in the rows
    pp_df[hydro_change, :el_MW] .= [
                                    3.0,    # kungfors
                                    91.0,   # olidans
                                    35.0,   # vargons
                                    46.0    # lilla edet
                                ]

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

function lines_props(
    lines_df::DataFrame
)

#=------------------------------------------------------------------------------
Function to add impendances and Slim to the raw OSM lines

values of R/km, X/km, and I_max in kA is currently based on
https://pandapower.readthedocs.io/en/latest/std_types/basic.html#lines
type 490-AL1/64-ST1A 110.0
susceptance is not modelled
use short line model
------------------------------------------------------------------------------=#

    # values assumed to be uniform
    r_per_km = 0.059 # Ω / km
    x_per_km = 0.37 # Ω / km
    z_per_km = r_per_km + im * x_per_km # Ω / kmB
    max_i_ka = 0.96 # kA
    Vnom = 130

    # resistance, reactance, impedance columns
    lines_df[!, :r_per_km] .= r_per_km
    lines_df[!, :x_per_km] .= x_per_km
    lines_df[!, :z_per_km] .= z_per_km

    # impedance due to parallel lines and length
    lines_df[!, :r_line] = [1 / sum(1/row[:r_per_km] for _ in 1:row[:circuits]) for row in eachrow(lines_df)]
    lines_df[!, :x_line] = [1 / sum(1/row[:x_per_km] for _ in 1:row[:circuits]) for row in eachrow(lines_df)]
    lines_df[!, :z_line] = [1 / sum(1/row[:z_per_km] for _ in 1:row[:circuits]) for row in eachrow(lines_df)]
    lines_df[!, :z_total] = lines_df[!, :z_line] .* lines_df[!, :length_km]

    # admittance of each lines
    # negative due to the convention in admittance matrix Y_ij equal to negative of admittance each line -y_ij
    lines_df[!, :y_total] = -1 ./ (lines_df[!, :z_total])

    # split into real imag parts for r,x,g,b
    lines_df[!, :r_total] = real(lines_df[!, :z_total])
    lines_df[!, :x_total] = imag(lines_df[!, :z_total])

    lines_df[!, :g_total] = real(lines_df[!, :y_total])
    lines_df[!, :b_total] = imag(lines_df[!, :y_total])

    # thermal limits
    lines_df[!, :max_i_ka] .= max_i_ka
    lines_df[!, :s_max] = sqrt(3) .* Vnom .* lines_df[!, :max_i_ka] .* lines_df[!, :circuits]

    # comparing to the admittance matrix method with incidence matrix
    Ybus, G, B = admittance_matrix(lines_df)

    return (; lines_df, Ybus, G, B)

end


function df_to_axisarrays(
    df::DataFrame
)

    # Convert the DataFrame to an Array
    matrix_df = Matrix(df)

    # Get column names and convert them to Symbols
    col_names = Symbol.(names(df))

    # Create an AxisArray with the column names as labels (for the second axis, i.e., columns)
    df_in_axis_arrays = AxisArrays.AxisArray(
                    matrix_df, 
                    1:size(matrix_df, 1), 
                    col_names
    )

    return df_in_axis_arrays

end


#=---------------------------------------------
POST - PROCESSING
---------------------------------------------=#

function sum_capacities(
    investment_df
)

    VGR_node = DataFrame(NODE = :VGR)
    total_inv = combine(investment_df, 
                            names(investment_df, Not(:NODE)) .=> sum .=> names(investment_df, Not(:NODE))
                            )

    total_inv = hcat(total_inv, VGR_node)
    investment_df = vcat(investment_df, total_inv)

    return investment_df

end


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


# Function to parse the .INC file
# convert GAMS EV .INC files into CSV
# this currently used for converting fleet availability .inc files
function parse_inc_file(inc_file::String)
    data = []
    open(inc_file, "r") do file
        for line in eachline(file)
            line = strip(line)
            # Skip comments and empty lines
            if isempty(line) || startswith(line, "*") # || startswith(line, "$")
                continue
            end
            # Assuming the data is space-separated or comma-separated
            row = split(line)
            if length(row) == 2
                push!(data, parse(Float64, row[2]))
            end
        end
    end
    return Float64.(data)
end

function combine_fleetava()

    # to adjust accordingly
    home_file = "fleetava_home_AGG.INC"  # Path to your .INC file
    data_home = parse_inc_file(home_file)

    h1_file = "fleetava_1h_AGG.INC"  # Path to your .INC file
    data_h1 = parse_inc_file(h1_file)

    h3_file = "fleetava_3h_AGG.INC"  # Path to your .INC file
    data_h3 = parse_inc_file(h3_file)

    h6_file = "fleetava_6h_AGG.INC"  # Path to your .INC file
    data_h6 = parse_inc_file(h6_file)

    fleet_df = DataFrame(
    :fleet_home => data_home,
    :fleet_h6 => data_h6,
    :fleet_h3 => data_h3,
    :fleet_h1 => data_h1,
    )

    CSV.write(joinpath(input_dir, "fleet_availability.csv"), fleet_df)

    return fleet_df

end

