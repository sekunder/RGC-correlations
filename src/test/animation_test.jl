include("../src/util/init.jl")
using JLD

D_real = load("/Users/sekunder/Documents/MATLAB/crcns_ret-1/cortex-backup/analysis/STRF/real/20080516_R1-1_STRF_127.jld")
D_sim = load("/Users/sekunder/Documents/MATLAB/crcns_ret-1/cortex-backup/analysis/STRF/sim/20080516_R1-1_simulated_127.jld")

strf = D_real["STRFs"][1]

(stat, msg) = animated_gif(strf, filename=joinpath(homedir(), "animation_test.gif"), fps=5, verbose=1)

#
STRFs_sim = load("/Users/sekunder/Documents/MATLAB/crcns_ret-1/cortex-backup/analysis/STRF/sim/20080516_R1-1_simulated_127.jld", "STRFs")
STRFs_real = load("/Users/sekunder/Documents/MATLAB/crcns_ret-1/cortex-backup/analysis/STRF/real/20080516_R1-1_STRF_127.jld", "STRFs")

for (idx, (s_r, s_s)) in enumerate(zip(STRFs_real, STRFs_sim))
    animated_gif(s_r, s_s; filename=joinpath(homedir(), "strf_comparison_$idx.gif"), fps=5, verbose=1, title=["Neuron $idx, real", "Neuron $idx, simulated"])
end
