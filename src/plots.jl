

Gen_Cap = CSV.read("Gen_Cap.csv", DataFrame)

Gen_Cap[:, :NODE]

VGR_gen = filter(:NODE => ==("VGR"), Gen_Cap)

Array(VGR[1, 2:end])

gen_names = names(VGR_gen[:, 2:end])

gen_caps = Array(VGR_gen[1, 2:end]) ./ 1000

f_gen, ax, bp = barplot(
                    gen_caps, 
                    color = 1:length(gen_caps), 
                    strokecolor = :black, 
                    strokewidth = 1,
                    width = 0.5,
                    gap = 0,
                    axis = (
                        xlabel="Generation Technologies",
                        ylabel="Capacity (GW)",
                        title="Generation Investment VGR",
                        limits=(nothing, nothing, 0, 14),
                        yticks=(0:2:14),
                        xticks=(1:length(gen_caps), gen_names),
                        xticklabelsize = 10,
                        xticklabelrotation = pi/4,
                        # xminorticks = 0.5:1:5.5,
                        # xminorticksvisible = true,
                        # xminorticksize = 8,
                        # topspinevisible=false,
                        # rightspinevisible=false,
                        # xgridvisible=false,
                        # ygridvisible=false,
                        # xticksvisible=false
                ),
                bar_labels = gen_caps,
                figure = (font = "Arial", size = (1024, 768))
            )
save("gen_inv.png", f_gen)

Sto_Cap = CSV.read("Sto_Inv.csv", DataFrame)

Sto_Cap[:, :NODE]

VGR_sto = filter(:NODE => ==("VGR"), Sto_Cap)

Array(VGR_sto[1, 2:end])

sto_names = names(VGR_sto[:, 2:end])

sto_caps = Array(VGR_sto[1, 2:end]) ./ 1000

f_sto, ax, bp = barplot(
                        sto_caps, 
                        color = 1:length(sto_caps), 
                        strokecolor = :black, 
                        strokewidth = 1,
                        width = 0.5,
                        gap = 0,
                        axis = (
                            xlabel="Storage Technologies",
                            ylabel="Capacity (GWh)",
                            title="Storage Investment VGR",
                            limits=(nothing, nothing, 0, nothing),
                            # yticks=(0:2:14),
                            xticks=(1:length(sto_caps), sto_names),
                            xticklabelsize = 10,
                            xticklabelrotation = pi/4,
                            # xminorticks = 0.5:1:5.5,
                            # xminorticksvisible = true,
                            # xminorticksize = 8,
                            # topspinevisible=false,
                            # rightspinevisible=false,
                            # xgridvisible=false,
                            # ygridvisible=false,
                            # xticksvisible=false
                    ),
                    bar_labels = sto_caps,
                    figure = (font = "Arial", size = (1024, 768))
                )
save("sto_inv.png", f_sto)
