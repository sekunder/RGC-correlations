
include("../util/init.jl")
# include("../util/metadata.jl")
include("../util/CRCNS/analysis_functions.jl")
using JLD, PyPlot

if !ispath(CRCNS_plots_dir)
    println("* Making path $CRCNS_plots_dir")
    mkpath(CRCNS_plots_dir)
end

now_time = now()
datestring = Dates.format(now_time, "YYYYmmdd")
timestring = Dates.format(now_time, "HHMMSS")

# plan: collect entropy values.
println("CRCNS_plot_I_2: Collecting entropy values that have been computed and plotting the resulting information ratios.")
I_2_real,I_2_simu = CRCNS_collect_entropy(;verbose=2)
# sorted_keys = sort(collect(keys(D)))
#
# I_2_real = Dict{Int,Vector{Float64}}()
# I_2_simu = Dict{Int,Vector{Float64}}()
# for k in keys(D)
#     I_2_real[k] = (D[k][1,:] - D[k][2,:]) ./ (D[k][1,:] - D[k][3,:])
#     I_2_simu[k] = (D[k][4,:] - D[k][5,:]) ./ (D[k][4,:] - D[k][6,:])
# end
sorted_keys = sort(collect(keys(I_2_real)))

println("* Creating figure")
fig = figure("I_2 ratio from available data", tight_layout=true, figsize=(10,10))
for k in keys(I_2_real)
    scatter(fill(k, length(I_2_real[k])), I_2_real[k], color="red", alpha=0.5, s=10, marker="o")
    scatter(fill(k, length(I_2_simu[k])), I_2_simu[k], color="blue", alpha=0.5, s=10, marker="o")
end
# plot(sort(keys(I_2_real)), [mean(I_2_real[k]) for k in sort(keys(I_2_real))], label="Real", color="red")
# plot(sort(keys(I_2_simu)), [mean(I_2_simu[k]) for k in sort(keys(I_2_simu))], label="Simulated", color="blue")
plot(sorted_keys, [mean(I_2_real[k]) for k in sorted_keys], label="mean Real", color="red")
plot(sorted_keys, [mean(I_2_simu[k]) for k in sorted_keys], label="mean Sim", color="blue")
title("Information ratio I_2")
xticks(sorted_keys, sorted_keys)
yticks(0.0:0.1:1.0, 0.0:0.1:1.0)
xlabel("Sample size")
ylabel(L"I_2")
grid("on")
legend()
println("  Saving figure to $(joinpath(CRCNS_plots_dir,"I_2_all-$datestring-$timestring.png"))")
savefig(joinpath(CRCNS_plots_dir,"I_2_all-$datestring-$timestring.png"))
close(fig)
