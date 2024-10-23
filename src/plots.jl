using DataFrames,CSV
import CairoMakie

# Plot Capacities in BarPlot
# Generation
Gen_Cap = CSV.read("Gen_Cap.csv", DataFrame)
Gen_Cap[:, :NODE]
VGR_gen = filter(:NODE => ==("VGR"), Gen_Cap)

gen_names = names(VGR_gen[:, 2:end])
gen_caps = Array(VGR_gen[1, 2:end])

f_gen, ax, bp = CairoMakie.barplot(
                    gen_caps, 
                    color = 1:length(gen_caps), 
                    strokecolor = :black, 
                    strokewidth = 1,
                    width = 0.5,
                    gap = 0,
                    axis = (
                        xlabel="Generation Technologies",
                        ylabel="Capacity (MW)",
                        title="Generation Investment VGR",
                        limits=(nothing, nothing, 0, nothing),
                        # yticks=(0:2:14),
                        xticks=(1:length(gen_caps), gen_names),
                        xticklabelsize = 10,
                        xticklabelrotation = pi/4,
                        ytickformat = "{:.0f}",
                        # xminorticks = 0.5:1:5.5,
                        # xminorticksvisible = true,
                        # xminorticksize = 8,
                        # topspinevisible=false,
                        # rightspinevisible=false,
                        # xgridvisible=false,
                        # ygridvisible=false,
                        # xticksvisible=false
                ),
                # bar_labels = gen_caps,
                figure = (font = "Arial", size = (600, 400))
            )
CairoMakie.save("gen_inv.png", f_gen)

# Storage Technologies
Sto_Cap = CSV.read("Sto_Inv.csv", DataFrame)
Sto_Cap[:, :NODE]
VGR_sto = filter(:NODE => ==("VGR"), Sto_Cap)

sto_names = names(VGR_sto[:, 2:end])
sto_caps = Array(VGR_sto[1, 2:end])

f_sto, ax, bp = CairoMakie.barplot(
                        sto_caps, 
                        color = 1:length(sto_caps), 
                        strokecolor = :black, 
                        strokewidth = 1,
                        width = 0.5,
                        gap = 0,
                        axis = (
                            xlabel="Storage Technologies",
                            ylabel="Capacity (MWh)",
                            title="Storage Investment VGR",
                            limits=(nothing, nothing, 0, nothing),
                            # yticks=(0:100000:300000),
                            xticks=(1:length(sto_caps), sto_names),
                            xticklabelsize = 10,
                            xticklabelrotation = pi/4,
                            ytickformat = "{:.0f}",
                            # xminorticks = 0.5:1:5.5,
                            # xminorticksvisible = true,
                            # xminorticksize = 8,
                            # topspinevisible=false,
                            # rightspinevisible=false,
                            # xgridvisible=false,
                            # ygridvisible=false,
                            # xticksvisible=false
                    ),
                    # bar_labels = sto_caps,
                    figure = (font = "Arial", size = (600, 400))
                )
CairoMakie.save("sto_inv.png", f_sto)


# Gen_Dispatch = CSV.read("Gen_Disp.csv", DataFrame)

# VGR_gen_dispatch = combine(groupby(Gen_Dispatch, [:Tech, :Period]), :Dispatch => sum => :Dispatch)
# VGR_gen_dispatch[!, :Node] .= "VGR"
# VGR_gen_dispatch = select(VGR_gen_dispatch, :Node, :Tech, :Period, :Dispatch)
# Gen_Dispatch = vcat(Gen_Dispatch, VGR_gen_dispatch)

# VGR_WOFF = filter([:Node, :Tech] => (x, y) -> x == "VGR" && y == "WOFF", Gen_Dispatch)
# VGR_AEC = filter([:Node, :Tech] => (x, y) -> x == "VGR" && y == "AEC", Gen_Dispatch)

# f_line, ax, line_plot = CairoMakie.lines(
#                                 1:504, # x
#                                 [VGR_WOFF[1:504, :Dispatch], 
#                                 VGR_PV[1:504, :Dispatch]],
#                                 axis = (
#                                     xlabel="Hour",
#                                     ylabel="Dispatch (MW)",
#                                     title="Generation of WOFF",
#                                     limits=(0, 500, 0, nothing),
#                                     # yticks=(0:2:14),
#                                     # xticks=(1:length(sto_caps), sto_names),
#                                     # xticklabelsize = 10,
#                                     # xticklabelrotation = pi/4,
#                                     # xminorticks = 0.5:1:5.5,
#                                     # xminorticksvisible = true,
#                                     # xminorticksize = 8,
#                                     # topspinevisible=false,
#                                     # rightspinevisible=false,
#                                     # xgridvisible=false,
#                                     # ygridvisible=false,
#                                     # xticksvisible=false
#                                 ),
#                                 figure = (font = "Arial", size = (600, 400))
# )

# Plot Dispatch in LinesPlot
# Generation
Gen_Dispatch = CSV.read("Gen_Disp.csv", DataFrame)

VGR_WOFF = filter([:Node, :Tech] => (x, y) -> x == "VGR" && y == "WOFF", Gen_Dispatch)
VGR_AEC = filter([:Node, :Tech] => (x, y) -> x == "VGR" && y == "AEC", Gen_Dispatch)

fig = CairoMakie.Figure(font = "Arial", size = (600, 400))
ax = CairoMakie.Axis(
    fig[1,1],
    xlabel="Hour",
    ylabel="Dispatch (MW)",
    title="Generation Dispatch",
    limits=(0, 500, 0, nothing),
    # yticks=(0:2:14),
    # xticks=(1:length(sto_caps), sto_names),
    # xticklabelsize = 10,
    # xticklabelrotation = pi/4,
    # xminorticks = 0.5:1:5.5,
    # xminorticksvisible = true,
    # xminorticksize = 8,
    # topspinevisible=false,
    # rightspinevisible=false,
    # xgridvisible=false,
    # ygridvisible=false,
    # xticksvisible=false
)
CairoMakie.lines!(ax, 1:504, VGR_WOFF[1:504, :Dispatch], label = "WOFF", color = :red)
CairoMakie.lines!(ax, 1:504, VGR_AEC[1:504, :Dispatch], label = "AEC", color = :blue)
CairoMakie.axislegend()
# fig[1, 2] = CairoMakie.Legend(fig, ax, "Technology", framevisible = true)
fig
CairoMakie.save("gen_disp.png", fig)

# Storage
Sto_Lv = CSV.read("Sto_Lv.csv", DataFrame)

VGR_LRC = filter([:Node, :Tech] => (x, y) -> x == "VGR" && y == "LRC", VGR_sto_lv)
VGR_BAT = filter([:Node, :Tech] => (x, y) -> x == "VGR" && y == "LI_EN", VGR_sto_lv)

fig = CairoMakie.Figure(font = "Arial", size = (600, 400))
ax = CairoMakie.Axis(
    fig[1,1],
    xlabel="Hour",
    ylabel="Storage Level (MWh)",
    title="Storage Level",
    limits=(0, 500, 0, nothing),
    ytickformat = "{:.0f}"
    # yticks=(0:2:14),
    # xticks=(1:length(sto_caps), sto_names),
    # xticklabelsize = 10,
    # xticklabelrotation = pi/4,
    # xminorticks = 0.5:1:5.5,
    # xminorticksvisible = true,
    # xminorticksize = 8,
    # topspinevisible=false,
    # rightspinevisible=false,
    # xgridvisible=false,
    # ygridvisible=false,
    # xticksvisible=false
)
CairoMakie.lines!(ax, 1:504, VGR_LRC[1:504, :Level], label = "LRC", color = :red)
CairoMakie.lines!(ax, 1:504, VGR_BAT[1:504, :Level], label = "Li-Ion", color = :blue)
CairoMakie.axislegend(position=:lt)
# fig[1, 2] = CairoMakie.Legend(fig, ax, "Technology", framevisible = true)
fig
CairoMakie.save("sto_lv.png", fig)


Q_Dispatch = CSV.read("Q_Disp.csv", DataFrame)

VGR_Q_AEC = filter([:Node, :Tech] => (x, y) -> x == "VGR" && y == "AEC", Q_Dispatch)
VGR_BAT = filter([:Node, :Tech] => (x, y) -> x == "VGR" && y == "LI_EN", Q_Dispatch)