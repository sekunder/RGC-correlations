include("../util/init.jl")
using JLD

# D = load("/Users/sekunder/Documents/MATLAB/crcns_ret-1/cortex-backup/analysis/information/20080516_R1-2.jld")
#
# D_real = load("/Users/sekunder/Documents/MATLAB/crcns_ret-1/cortex-backup/analysis/STRF/real/20080516_R2-3_STRF_17179869183.jld")
D_sim = load("/Users/sekunder/Documents/MATLAB/crcns_ret-1/cortex-backup/analysis/STRF/sim/20080516_R2-3_simulated_17179869183.jld")
spikes_sim = load(joinpath(CRCNS_STRF_dir,"sim","20080516_R2-3_simulated_17179869183.jld"),"spikes")

X_sim_33554431 = transpose(raster(spikes_sim, 0.020))

P_2_33554431 = second_order_model(X_sim_33554431, 1:25; verbose=2)
