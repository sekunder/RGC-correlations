help_string = """
$(basename(@__FILE__))

Computes STRFs for CRCNS data, then generates simulated data by presenting
stimulus to those STRFs.

Usage:
  [julia] $(basename(@__FILE__)) --defaults | --dir data_dir | file1 file2 ... [options]
  Processes files specified by first argument/flags:
    --defaults : read all .mat files in CRCNS_Data_dir, also sets output directories
    --dir data_dir : real all .mat files in dir
    file1 file2 ... : read files in list from CRCNS_Data_dir
  Options:
    --verbose 0|1|2 : controls verbosity of output
    --output_dir_real path, --output_dir_sim path : where to write .jld files
    --skip file1 file2 ... : which files to skip (including overriding defaults)
    --poisson_time t, --poisson_loops n, --poisson_spikes m : Stop poisson process generation if it lasts longer than t seconds, n loops, or m spikes.
        By default, uses 600.0 s (no value for loops, spikes)
    --help : display this message and exit
"""
# CRCNS_script_version=v"0.3"
# plan: for each mat file, for each recording index:
# 1. load the stimulus, call CRCNS_output_STRFs to save them to file
# 2. for each STRF: compute response to stimulus, scale response, generate spike trains using scaled response.
include("../util/init.jl")
using MAT, JLD
# include("../probability/probability.jl")
# include("../spikes/spikes.jl")
# include("../stimulus/stimulus.jl")
include("../util/CRCNS/analysis_functions.jl")
include("../util/nonlinearities.jl")
include("../spikes/CRCNS/CRCNS.jl")
include("../stimulus/CRCNS/CRCNS_Stimulus.jl")
# include("../util/misc.jl")
# include("../util/constants.jl")
# using Spikes, Stimulus, Probability

cline_args = process_args(ARGS; parse_flags=["verbose", "poisson_time","poisson_loops","poisson_spikes"], bool_flags=["help","defaults"])
verbose = get(cline_args,"verbose",[0])[1]
poisson_time = Float64(get(cline_args,"poisson_time",[600.0])[1])
poisson_loops = get(cline_args,"poisson_loops",[0])[1]
poisson_spikes = get(cline_args,"poisson_spikes",[0])[1]

if cline_args["help"]
    println(help_string)
    println("VERSION: $CRCNS_script_version")
    exit(0)
end

println("-" ^ 80)
println("CRCNS_generate_STRFs $CRCNS_script_version")
# println(help_string)

data_dir = get(cline_args,"dir",[CRCNS_data_dir])[1]
output_dir_real = get(cline_args, "output_dir_real", [joinpath(CRCNS_STRF_dir, "real")])[1]
output_dir_sim = get(cline_args, "output_dir_sim", [joinpath(CRCNS_STRF_dir, "sim")])[1]
file_list = haskey(cline_args,"dir") ? readdir(cline_args["dir"][1]) : cline_args["0"]

if cline_args["defaults"]
    data_dir = joinpath(CRCNS_dir, "Data")
    file_list = readdir(data_dir)
    output_dir_real = joinpath(CRCNS_STRF_dir, "real")
    output_dir_sim = joinpath(CRCNS_STRF_dir, "sim")
end

file_list = setdiff(file_list, get(cline_args,"skip",String[]))

for dir in [data_dir, output_dir_real, output_dir_sim]
    if !isdir(dir)
        println("* Creating path $dir")
        mkpath(dir)
    end
end

mat_files = filter(x -> endswith(x, ".mat"), file_list)
real_jld_files = filter(x -> endswith(x, ".jld"), readdir(output_dir_real))
sim_jld_files = filter(x -> endswith(x, ".jld"), readdir(output_dir_sim))

println("Verbosity level: $verbose")
println("Poisson options:")
println("\tReal time: $poisson_time")
println("\tLoops    : $poisson_loops")
println("\tSpikes   : $poisson_spikes")
println("Preparing to process these files:")
println("\t$(join(mat_files,"\n\t"))")

@everywhere function process_file(mat_file)
    lf = open(replace(abspath(@__FILE__), ".jl", ".$(myid()).log"), "a")
    println(lf, "* Processing file $mat_file")
    # poor planning on my part means I still have to open the damn files here to
    # get the number of recording indexes. Oh well.
    vars = try
            matread(joinpath(data_dir, mat_file))
        catch y
            println(lf, "! Error occurred, skipping file.")
            # show(y)
            println(lf, y)
            rethrow(y)
        end
    recordings = 1:length(vars["datainfo"]["RecNo"])
    # recordings = [2]
    vars = 0 # just a hint for garbage collecion
    println(lf, "  Found $(length(recordings)) recordings")
    for rec_idx in recordings
        println(lf, "  Recording $rec_idx")
        # real_files = filter(x -> startswith(x, "$(remove_extension(mat_file))-$rec_idx"), real_jld_files)
        # sim_files = filter(x -> startswith(x, "$(remove_extension(mat_file))-$rec_idx"), sim_jld_files)
        # for (rf, sf) in zip(real_files,sim_files)
        #     cv_r = read(joinpath(output_dir_real, rf), "CRCNS_script_version")
        #     if VersionNumber(cv_r) < CRCNS_script_version
        #         println("  Found file $rf using script version $cv_r. Skipping processing.")
        #         continue
        #     end
        #     cv_s = read(joinpath(output_dir_sim, sf), "CRCNS_script_version")
        #     if VersionNumber(cv_s) < CRCNS_script_version
        #         println("  Found file $sf using script version $cv_s. Skipping processing.")
        #         continue
        #     end
        # end
        print(lf, "    Computing real STRFs...")

        stim, spikes, spike_hist, STRFs = CRCNS_output_STRFs(joinpath(data_dir, mat_file), rec_idx, output_dir_real;
            CRCNS_script_version=CRCNS_script_version, verbose=verbose, single_rec=(length(recordings) == 1))
        println(lf, "done")
        L = zeros(spike_hist)
        ST_simulated = Vector{Vector{Float64}}(n_cells(spikes))
        ST_status = Vector{Vector{Symbol}}(n_cells(spikes))
        println(lf, "    Simulating $(n_cells(spikes)) cells, using phi = sigmoid, Q = norm(r * tau - n) / N_frames")
        println(lf, "    [r = response computed | s = response scaled to match STRFs | p = poisson process spike train generated]")
        print(lf, "      ")
        for (idx, RF) in enumerate(STRFs)
            print(lf, "$idx[")
            (r,tau) = STRF_response(RF, stim, flip_STRF_time=true)
            n = spike_hist[:,idx]
            print(lf, "r|")

            c_range = (1.0:0.1:maximum(n)) / tau
            h_range = [10.0^k for k in 2.0:0.1:4.0]
            x0_range = decimal_round(minimum(r),2):0.001:decimal_round(maximum(r),2)
            theta_ranges = [c_range, h_range, x0_range]
            L[:,idx], theta_opt, Q_opt = scale_response(r, n, sigmoid, (u,v) -> (norm(u * tau - v) / length(u)); d=3, ranges=theta_ranges, save_fun=Float64[])
            print(lf, "s|")

            if verbose > 1
                expected_spikes = sum_kbn((L[:,idx] * tau)[:])
                print(lf, "[Expected #spikes = $expected_spikes]")
            end

            ST_status[idx] = Vector{Symbol}()
            time_arr = zeros(1)
            ST_simulated[idx] = inhomogeneous_poisson_process(L[:,idx], tau;
                sampling_rate_factor=10, max_real_time=poisson_time, max_loops=poisson_loops, max_spikes=poisson_spikes,
                exit_status=ST_status[idx], total_time=time_arr)
            print("p")
            if verbose > 1
                print(lf, "[Got $(length(ST_simulated[idx])) spikes in $(time_arr[1]) s]")
            end
            print(lf, "] ")
        end

        print(lf, "\n    Computing simulated STRFs...")
        sim_spikes = SpikeTrains(ST_simulated, spikes.I;
            comment="Simulated spike train for CRCNS data $mat_file, recording index $rec_idx",
            CRCNS_script_version=CRCNS_script_version,
            poisson_status=ST_status, poisson_rate=L)
        hide_metadata!(sim_spikes, :poisson_status)
        hide_metadata!(sim_spikes, :poisson_rate)
        sim_hist = histogram(sim_spikes, frame_time(stim); N_bins=n_frames(stim))
        sim_STRFs = compute_STRFs(sim_hist, stim)
        println(lf, "done")
        if verbose > 1
            println(lf, sim_spikes)
        end

        indexes = index_set_to_int(sim_spikes.I)
        sim_filename = "$(remove_extension(mat_file))-$(rec_idx)_simulated_$indexes.jld"
        println(lf, "    Writing simulated spike trains and computed STRFs to $(joinpath(output_dir_sim,sim_filename))")
        # save(joinpath(output_dir_sim, sim_filename), "CRCNS_script_version", CRCNS_script_version, "timestamp", now(), "STRFs", sim_STRFs, "spikes", sim_spikes)
        jldopen(joinpath(output_dir_sim, sim_filename), "w") do file
            addrequire(file, Spikes)
            addrequire(file, GrayScaleStimuli)
            write(file, "CRCNS_script_version", CRCNS_script_version)
            write(file, "timestamp", now())
            write(file, "STRFs", sim_STRFs)
            write(file, "spikes", sim_spikes)
        end
    end
    close(lf)
end

pmap(process_file, mat_files)

# for mat_file in mat_files
#     println("* Processing file $mat_file")
#     # poor planning on my part means I still have to open the damn files here to
#     # get the number of recording indexes. Oh well.
#     vars = try
#             matread(joinpath(data_dir, mat_file))
#         catch y
#             println("! Error occurred, skipping file.")
#             show(y)
#             continue
#         end
#     recordings = 1:length(vars["datainfo"]["RecNo"])
#     # recordings = [2]
#     vars = 0 # just a hint for garbage collecion
#     println("  Found $(length(recordings)) recordings")
#     for rec_idx in recordings
#         println("  Recording $rec_idx")
#         real_files = filter(x -> startswith(x, "$(remove_extension(mat_file))-$rec_idx"), real_jld_files)
#         sim_files = filter(x -> startswith(x, "$(remove_extension(mat_file))-$rec_idx"), sim_jld_files)
#         for (rf, sf) in zip(real_files,sim_files)
#             cv_r = read(joinpath(output_dir_real, rf), "CRCNS_script_version")
#             if VersionNumber(cv_r) < CRCNS_script_version
#                 println("  Found file $rf using script version $cv_r. Skipping processing.")
#                 continue
#             end
#             cv_s = read(joinpath(output_dir_sim, sf), "CRCNS_script_version")
#             if VersionNumber(cv_s) < CRCNS_script_version
#                 println("  Found file $sf using script version $cv_s. Skipping processing.")
#                 continue
#             end
#         end
#         # if !isempty(real_files) && !isempty(sim_files)
#         #     println("  Found files")
#         #     println("    Real: $(join(real_files,","))")
#         #     println("    Sim : $(join(sim_files,","))")
#         #     println("  Skipping processing.")
#         #     continue
#         # end
#         print("    Computing real STRFs...")
#
#         stim, spikes, spike_hist, STRFs = CRCNS_output_STRFs(joinpath(data_dir, mat_file), rec_idx, output_dir_real;
#             CRCNS_script_version=CRCNS_script_version, verbose=verbose, single_rec=(length(recordings) == 1))
#         println("done")
#         L = zeros(spike_hist)
#         ST_simulated = Vector{Vector{Float64}}(n_cells(spikes))
#         ST_status = Vector{Vector{Symbol}}(n_cells(spikes))
#         println("    Simulating $(n_cells(spikes)) cells, using phi = sigmoid, Q = norm(r * tau - n) / N_frames")
#         println("    [r = response computed | s = response scaled to match STRFs | p = poisson process spike train generated]")
#         print("      ")
#         for (idx, RF) in enumerate(STRFs)
#             print("$idx[")
#             (r,tau) = STRF_response(RF, stim, flip_STRF_time=true)
#             n = spike_hist[:,idx]
#             print("r|")
#
#             c_range = (1.0:0.1:maximum(n)) / tau
#             h_range = [10.0^k for k in 2.0:0.1:4.0]
#             x0_range = decimal_round(minimum(r),2):0.001:decimal_round(maximum(r),2)
#             theta_ranges = [c_range, h_range, x0_range]
#             L[:,idx], theta_opt, Q_opt = scale_response(r, n, sigmoid, (u,v) -> (norm(u * tau - v) / length(u)); d=3, ranges=theta_ranges, save_fun=Float64[])
#             print("s|")
#
#             if verbose > 1
#                 expected_spikes = sum_kbn((L[:,idx] * tau)[:])
#                 print("[Expected #spikes = $expected_spikes]")
#             end
#
#             ST_status[idx] = Vector{Symbol}()
#             time_arr = zeros(1)
#             ST_simulated[idx] = inhomogeneous_poisson_process(L[:,idx], tau;
#                 sampling_rate_factor=10, max_real_time=poisson_time, max_loops=poisson_loops, max_spikes=poisson_spikes,
#                 exit_status=ST_status[idx], total_time=time_arr)
#             print("p")
#             if verbose > 1
#                 print("[Got $(length(ST_simulated[idx])) spikes in $(time_arr[1]) s]")
#             end
#             print("] ")
#         end
#
#         print("\n    Computing simulated STRFs...")
#         sim_spikes = SpikeTrains(ST_simulated, spikes.I;
#             comment="Simulated spike train for CRCNS data $mat_file, recording index $rec_idx",
#             CRCNS_script_version=CRCNS_script_version,
#             poisson_status=ST_status, poisson_rate=L)
#         hide_metadata!(sim_spikes, :poisson_status)
#         hide_metadata!(sim_spikes, :poisson_rate)
#         sim_hist = histogram(sim_spikes, frame_time(stim); N_bins=n_frames(stim))
#         sim_STRFs = compute_STRFs(sim_hist, stim)
#         println("done")
#         if verbose > 1
#             println(sim_spikes)
#         end
#
#         indexes = index_set_to_int(sim_spikes.I)
#         sim_filename = "$(remove_extension(mat_file))-$(rec_idx)_simulated_$indexes.jld"
#         println("    Writing simulated spike trains and computed STRFs to $(joinpath(output_dir_sim,sim_filename))")
#         # save(joinpath(output_dir_sim, sim_filename), "CRCNS_script_version", CRCNS_script_version, "timestamp", now(), "STRFs", sim_STRFs, "spikes", sim_spikes)
#         jldopen(joinpath(output_dir_sim, sim_filename), "w") do file
#             addrequire(file, Spikes)
#             addrequire(file, GrayScaleStimuli)
#             write(file, "CRCNS_script_version", CRCNS_script_version)
#             write(file, "timestamp", now())
#             write(file, "STRFs", sim_STRFs)
#             write(file, "spikes", sim_spikes)
#         end
#     end
# end

println("CRCNS_generate_STRFs: END SCRIPT $(now())")
println("-" ^ 80)
