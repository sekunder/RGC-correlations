help_string = """
$(@__FILE__)

Computes STRFs for CRCNS data, taking advantage of multiple available cores.

Usage:
    [julia] $(@__FILE__) --defaults | --dir data_dir | file1 file2 ... [options]
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

@everywhere include("../util/init.jl")
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
println("CRCNS_generate_STRFs $CRCNS_script_version")
println("For help, run")
println("    [julia] $(@__FILE__) --help")

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
