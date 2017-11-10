
CRCNS_script_version=v"0.1"
# plan: for each mat file, for each recording index:
# 1. load the stimulus, call CRCNS_output_STRFs to save them to file
# 2. for each STRF: compute response to stimulus, scale response, generate spike trains using scaled response.

using MAT, JLD
include("../probability/probability.jl")
include("../spikes/spikes.jl")
include("../stimulus/stimulus.jl")
include("../util/CRCNS/analysis_functions.jl")
include("../util/nonlinearities.jl")
include("../util/misc.jl")
using Spikes, Stimulus, Probability

data_dir = joinpath(Stimulus.default_CRCNS_dir, "Data")
output_dir_real = joinpath(Stimulus.default_CRCNS_dir, "analysis", "STRF", "real")
output_dir_sim = joinpath(Stimulus.default_CRCNS_dir, "analysis", "STRF", "sim")

mat_files = filter(x -> endswith(x, ".mat"), readdir(data_dir))
# mat_files = filter(x -> x == "20080516_R1.mat", readdir(data_dir)) #for debugging

println("-" ^ 80)
println("CRCNS_generate_STRFs: BEGIN SCRIPT")
for dir in [data_dir, output_dir_real, output_dir_sim]
    if !isdir(dir)
        println("* Creating path $dir")
        mkpath(dir)
    end
end
println("Computing STRFs for CRCNS data, then generating simulated data by presenting stimulus to those STRFs.")
println("Preparing to process these files:")
println(mat_files)

for mat_file in mat_files
    println("* Processing file $mat_file")
    # poor planning on my part means I still have to open the damn files here to
    # get the number of recording indexes. Oh well.
    vars = matread(joinpath(data_dir, mat_file))
    recordings = 1:length(vars["datainfo"]["RecNo"])
    # recordings = [2]
    vars = 0 # just a hint for garbage collecion
    for rec_idx in recordings
        println("  Recording $rec_idx")
        print("    Computing real STRFs...")
        # stim = CRCNS_Stimulus(joinpath(Stimulus.default_CRCNS_dir, data_dir, mat_file), rec_idx)
        # spike_hist = histogram(spikes, frame)
        stim, spikes, spike_hist, STRFs = CRCNS_output_STRFs(joinpath(data_dir, mat_file), rec_idx, output_dir_real; CRCNS_script_version=CRCNS_script_version)
        println("done")
        L = zeros(spike_hist)
        ST_simulated = Vector{Vector{Float64}}(n_cells(spikes))
        println("    Simulating, using phi = sigmoid, Q = norm(r * tau - n) / N_frames")
        print("      ")
        for (idx, RF) in enumerate(STRFs)
            print("$idx[")
            (r,tau) = STRF_response(RF, stim, flip_STRF_time=true)
            n = spike_hist[:,idx]
            print("response|")

            c_range = (1.0:0.1:maximum(n)) / tau
            h_range = [10.0^k for k in 2.0:0.1:4.0]
            x0_range = decimal_round(minimum(r),2):0.001:decimal_round(maximum(r),2)
            theta_ranges = [c_range, h_range, x0_range]
            L[:,idx], theta_opt, Q_opt = scale_response(r, n, sigmoid, (u,v) -> (norm(u * tau - v) / length(u)); d=3, ranges=theta_ranges, save_fun=Float64[])
            print("scaled|")

            ST_simulated[idx] = inhomogeneous_poisson_process(L[:,idx], tau; sampling_rate_factor=10)
            print("poisson] ")
        end

        print("\n    Computing simulated STRFs...")
        sim_spikes = SpikeTrains(ST_simulated, spikes.I; comment="Simulated spike train for CRCNS data $mat_file, recording index $rec_idx")
        sim_hist = histogram(sim_spikes, frame_time(stim))
        sim_STRFs = compute_STRFs(sim_hist, stim)
        println("done")

        indexes = mapreduce(x -> 2^(x-1), +, sim_spikes.I)
        sim_filename = "$(remove_extension(mat_file))-$(rec_idx)_simulated_$indexes.jld"
        println("    Writing simulated spike trains and computed STRFs to $(joinpath(output_dir_sim,sim_filename))")
        save(joinpath(output_dir_sim, sim_filename), "CRCNS_script_version", CRCNS_script_version, "timestamp", now(), "STRFs", sim_STRFs, "spikes", sim_spikes)
    end
end

println("CRCNS_generate_STRFs: END SCRIPT")
println("-" ^ 80)
