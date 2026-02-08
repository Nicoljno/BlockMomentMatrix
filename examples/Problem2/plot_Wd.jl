using LaTeXStrings
using Plots
using DelimitedFiles


classic_arr_2=readdlm("plot_classical_W2.txt")
classic_arr_3=readdlm("plot_classical_W3.txt")
classic_arr_4=readdlm("plot_classical_W4.txt")
classic_arr_5=readdlm("plot_classical_W5.txt")
quantum_SIC_d3=readdlm("SIC_d3.txt")
classic_arr_x2 = classic_arr_2[:,1]/2
classic_arr_y2 = classic_arr_2[:,2] 
classic_arr_x3 = classic_arr_3[:,1]/2
classic_arr_y3 = classic_arr_3[:,2] 
classic_arr_x4 = classic_arr_4[:,1]/2
classic_arr_y4 = classic_arr_4[:,2] 
classic_arr_x5 = classic_arr_5[:,1]/2
classic_arr_y5 = classic_arr_5[:,2] 
SIC_x3 = quantum_SIC_d3[:,1]
SIC_y3 = quantum_SIC_d3[:,2] 

gr(size = (1000, 500), linewidth = 4)


my_palette = [
    colorant"#1f77b4",
    colorant"#8c6d31",
    colorant"#2ca25f",
    colorant"#B0290E",
    colorant"#5B10B5",
]

p = plot(
    classic_arr_y3,  
    classic_arr_x3,  
    palette=my_palette,
    legend = :bottomright,
    label = "",
    framestyle = :box,
    xlims    = (0, 0.08),
    ylims    = (0.71, 1.005),
    xticks = 0:0.02:0.08,
    left_margin=5.5Plots.mm,
    right_margin=5.8Plots.mm,
    bottom_margin=5Plots.mm,
    xlabel = L"ω", 
    ylabel = L"W_{\mathrm{d}}",
    grid = false,
    lw = 3,
    xtickfont=font(20), 
    ytickfont=font(20), 
    guidefont=font(26), 
    legendfont=font(20),
)

plot!(classic_arr_y4, classic_arr_x4, lw=3, label="")
plot!(classic_arr_y5, classic_arr_x5, lw=3, label="")
plot!(classic_arr_y2, classic_arr_x2, lw=3, label="")

plot!(twinx(p), SIC_y3, SIC_x3; lw=3, label="", ylabel = L"P_{\mathrm{s}}", xlims = (0, 0.08), ylims    = (0.32, 0.75), color=colorant"#5B10B5", 
    yforeground_color_axis = colorant"#5B10B5",
    yforeground_color_text = colorant"#5B10B5",
    yforeground_color_guide = colorant"#5B10B5",
    ytickfontcolor = colorant"#5B10B5",
    yforeground_color_border = colorant"#5B10B5",
    yguidefontcolor = colorant"#5B10B5",
    linestyle = :dash,
    ytickfont=font(20), 
    guidefont=font(26), 
    legendfont=font(20),
    )

ptop = twiny(p)
ticks = 0:0.02:0.08
plot!(ptop;
    xlims  = (0, 0.08),
    xticks = (collect(ticks), fill("", length(ticks))),
    xlabel = "",
    yaxis  = false,
    grid   = false,
    legend = false,
)
#hline!(p, [ylims(p)[2]]; lc=:black, lw=3, linestyle = :solid, label="")
hline!([1], linestyle=:dash, color=:gray, label=false)

xproxy = [1.0, 1.1]   # outside (0, 0.08)

plot!(xproxy, [0.0, 0.0],
      color = my_palette[4], lw = 1.5, linestyle = :solid,
      label = L"W_2")

plot!(xproxy, [0.0, 0.0],
      color = my_palette[1], lw = 1.5, linestyle = :solid,
      label = L"W_3")

plot!(xproxy, [0.0, 0.0],
      color = my_palette[2], lw = 1.5, linestyle = :solid,
      label = L"W_4")

plot!(xproxy, [0.0, 0.0],
      color = my_palette[3], lw = 1.5, linestyle = :solid,
      label = L"W_5")

plot!(xproxy, [0.0, 0.0],
      color = colorant"#5B10B5", lw = 1.0, linestyle = :dash,
      label = L"\textbf{E}_{\mathrm{Hesse}}")

savefig("plot_Wd.pdf")

