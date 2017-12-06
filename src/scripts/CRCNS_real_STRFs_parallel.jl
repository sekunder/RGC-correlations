help_string = """
$(@__FILE__)

Computes STRFs for CRCNS data, taking advantage of multiple available cores.

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
        --poisson_time t, --poisson_loops n, --poisson_spikes m : Stop poisson process generation if it lasts longer than t seconds, n loops, or m spikes
        --help : display this message and exit
"""

include("../util/init.jl")
@everywhere using MAT, JLD
@everywhere include("../util/CRCNS/analysis_functions.jl")
@everywhere include("../util/nonlinearities.jl")

cline_args = process_args(ARGS; parse_flags=["verbose", "poisson_time","poisson_loops","poisson_spikes"], bool_flags=["help","defaults"])
verbose = get(cline_args,"verbose",[0])[1]
poisson_time = Float64(get(cline_args,"poisson_time",[0.0])[1])
poisson_loops = get(cline_args,"poisson_loops",[0])[1]
poisson_spikes = get(cline_args,"poisson_spikes",[0])[1]

if cline_args["help"]
    println(help_string)
    println("VERSION: $CRCNS_script_version")
    exit(0)
end

println("-" ^ 80)
println("$(basename(@__FILE__)) $CRCNS_script_version")
println("For help, run")
println("    [julia] $(basename(@__FILE__)) --help")

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

# game plan: write functions that perform the following whole procedures:
# 1. (mat filename, rec_num) -> write STRF to JLD file
# 2. (jld filename, cell num) -> stimulate STRF, fit parameters, generate poisson spikes

# first, let's create a list of (mat filename, rec num) tuples for each mat file this script is supposed to process.
multi_rec_tuples = []
single_rec_tuples = []
println("* Gathering information on recordings")
for mat_file in mat_files
    vars = try
            matread(joinpath(data_dir, mat_file))
        catch y
            Dict("datainfo"=>Dict("RecNo"=>[]))
        end
    Nrec = length(vars["datainfo"]["RecNo"])
    if Nrec == 1
        push!(single_rec_tuples, (mat_file,1))
    elseif Nrec > 1
        for i = 1:Nrec
            push!(multi_rec_tuples,i)
        end
    end
    println("  Found $Nrec recordings in $mat_file")
end

@everywhere function process_multi_rec(mat_file, rec_idx)
    return CRCNS_output_STRFs(joinpath(data_dir, mat_file), rec_idx, output_dir_real;
            CRCNS_script_version=CRCNS_script_version, verbose=verbose, single_rec=false)
end

@everywhere function process_single_rec(mat_file, rec_idx)
    return CRCNS_output_STRFs(joinpath(data_dir, mat_file), rec_idx, output_dir_real;
            CRCNS_script_version=CRCNS_script_version, verbose=verbose, single_rec=true)
end

println("* Attempting to generate STRFs")

STRF_stuff_multi = fetch(pmap(process_multi_rec, multi_rec_tuples; on_error=ex->(print("! error encountered: "); show(ex))))
STRF_stuff_single = fetch(pmap(process_single_rec, single_rec_tuples; on_error=ex->(print("! error encountered: "); show(ex))))

println("* STRFs generated")

@everywhere function simulated_response(stim, spikes, spike_hist, STRFs)
    L = zeros(spike_hist)
    resp(idx, RF) = begin
            (r,tau) = STRF_response(RF, stim, flip_STRF_time = true)
            n = spike_hist[:,idx]

            c_range = (1.0:0.1:maximum(n)) / tau
            h_range = [10.0^k for k in 2.0:0.1:4.0]
            x0_range = decimal_round(minimum(r),2):0.001:decimal_round(maximum(r),2)
            theta_ranges = [c_range, h_range, x0_range]

            L, theta_opt, Q_opt = scale_response(r, n, sigmoid, (u,v) -> norm(u*tau - v) / length(u); d=3, ranges=theta_ranges, save_fun=Float64[])
            ST_status = Vector{Symbol}()
            time_arr=zeros(1)
            ST_simulated = inhomogeneous_poisson_process(L, tau;
                sampling_rate_factor=10, max_real_time=poisson_time, max_loops=poisson_loops, max_spikes=poisson_spikes,
                exit_status=ST_status[idx], total_time=time_arr)
            return L, ST_simulated, ST_status
        end
    return pmap(resp, enumerate(STRFs); on_error=ex->(println("! error with resp: "); show(ex)))
end

# @everywhere function simulated_response(idx, RF)
#     (r,tau) = STRF_response(RF, STRF_stuff[1], flip_STRF_time=true)
#     n = STRF_stuff[3][:,idx]
#     c_range = (1.0:0.1:maximum(n)) / tau
#     h_range = [10.0^k for k in 2.0:0.1:4.0]
#     x0_range = decimal_round(minimum(r),2):0.001:decimal_round(maximum(r),2)
#     theta_ranges = [c_range, h_range, x0_range]
#     L, theta_opt, Q_opt = scale_response(r, n, sigmoid, (u,v)->norm(u*tau - v)/length(u); d=3, ranges=theta_ranges, save_fun=Float64[])
#
# end
