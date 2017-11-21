help_string = """
$(@__FILE__)

Fits probability models to subsamples of real and simulated data, then writes
these models to .jld files.

Usage:
  [julia] $(@__FILE__) --defaults | --dir sim_jld_dir | file1 file2 ... [options]
  Processes files specified by opening flags/files list:
    --defaults : Attempts to process every .jld file in the "sim" directory under `CRCNS_STRF_dir`
    --dir sim_jld_dir : Attempts to process every .jld file in `sim_jld_dir`
    file1 file2 ... : Attemps to process each specified file, looking in `CRCNS_STRF_dir`/sim
  Options:
    --verbose 0|1|2 : controls verbosity of output.
    --n_trials n : Run `n` trials for each sample size. Default 20
    --bin_size t : Compute spike rasters with bin size `t`. If `t` parses as an integer, assumes `t` is in milliseconds. If `t` is parsed as a float it is time in seconds.
    --maxtime t : Maxtime option for NLopt
    --data_dir dir : Specifies location of .mat files. Default is `CRCNS_data_dir`
    --output_dir dir : specifies the output directory. Default is `CRCNS_information_dir`
    --help : display this message and exit
"""

include("../util/init.jl")
include("../util/CRCNS/analysis_functions.jl")
using JLD

# plan: loop through folder of simulated STRFs (those files also have spike
# trains!):

# 1. load the simulated and real spike trains. compute rasters.

# 2. loop through: sample sizes, number of trials. Fit P_*, compute entropy (if
# appropriate), add this to a dictionary structure, etc.

# cline_args = process_args(ARGS, parse_defaults=Dict("n_trials"=>20, "bin_size"=>10e-3, "maxtime"=>0))
cline_args = process_args(ARGS, parse_flags=["n_trials", "bin_size","maxtime"], bool_flags=["help","defaults"])
verbose = get(cline_args,"verbose",[0])[1]
n_trials = get(cline_args,"n_trials",[20])[1]
# To make life easy, I want to be able to say julia /script 10 to do 10ms
temp_bin_size = get(cline_args,"bin_size",[10e-3])[1]
bin_size = isa(temp_bin_size, Integer) ? temp_bin_size / 1000.0 : temp_bin_size
maxtime = get(cline_args,"maxtime",[0])[1]

default_size_range = 5:5:40

if cline_args["help"]
    println(help_string)
    println("VERSION: $CRCNS_script_version")
    exit(0)
end

println("-" ^ 80)
println("CRCNS_output_information $CRCNS_script_version: BEGIN SCRIPT $(now())")

sim_jld_dir = get(cline_args,"dir",[joinpath(CRCNS_STRF_dir, "sim")])[1]
information_dir = get(cline_args,"output_dir",[CRCNS_information_dir])[1]
data_dir = get(cline_args,"data_dir",[CRCNS_data_dir])[1]

for dir in [sim_jld_dir,CRCNS_information_dir,data_dir]
    if !isdir(dir)
        println("* Creating path $dir")
        mkpath(dir)
    end
end

file_list = haskey(cline_args,"dir") ? readdir(cline_args["dir"][1]) : cline_args["0"]

if cline_args["defaults"]
    sim_jld_dir = joinpath(CRCNS_STRF_dir,"sim")
    file_list = readdir(sim_jld_dir)
    information_dir = CRCNS_information_dir
    data_dir = CRCNS_data_dir
end

sim_jld_files = filter(x -> endswith(x,".jld"), file_list)

println("Loading simulated and real spikes, then fitting P_1, P_2 to subsamples of the data, using bin size $bin_size s")
if maxtime > 0
    println("NLopt maxtime = $maxtime")
end
println("Will run $n_trials trials for each sample size")
println("Default size range is $default_size_range")
println("Will write output to $information_dir")

println("Preparing to process these files:")
println("\t$(join(sim_jld_files,"\n\t"))")

successful_files = String[]

for sim_file in sim_jld_files
    println("* Processing file $sim_file")
    root_name = sim_file[1:11] # grab the "2008XXXX_RX" portion of the filename
    rec_idx = 0
    try
        # grab the recording index.
        rec_idx = parse(Int, sim_file[13])
    catch y
        if isa(y, ArgumentError)
            rec_idx = 0
        else
            rethrow(y)
        end
    end
    if rec_idx == 0
        println("! Found $sim_file, but could not determine recording index. Skipping file")
        continue
    end
    if !isfile(joinpath(data_dir, "$root_name.mat"))
        println("! Found $sim_file, but could not find $root_name.mat. Skipping file.")
        continue
    end
    println("* Processing $root_name, recording index $rec_idx")
    println("  Loading CRCNS spikes from $root_name.mat")
    # what the fuck is up with how julia handles global variables?! without the
    # line below, this script gives me an UndefVarError for real_spikes, even
    # though the verbose output shows it successfully ran. What's worse is this
    # will probably cause some overhead because it's changing the type so
    # drastically.
    real_spikes = 0
    sim_spikes = 0
    try
        real_spikes = Spikes.CRCNS_get_spikes_from_file(joinpath(data_dir, "$root_name.mat"), rec_idx)
    catch y
        println("! Exception occurred. Skipping file.")
        show(y)
        continue
    end
    println("  Loading simulated spikes from $sim_file")
    try
        sim_spikes = load(joinpath(sim_jld_dir, sim_file), "spikes")
    catch y
        println("! Exception occurred. Skipping file.")
        show(y)
        continue
    end
    if n_cells(real_spikes) != n_cells(sim_spikes)
        println("! CRCNS spikes for $(n_cells(real_spikes)) cells, but simulated data for $(n_cells(sim_spikes)). Skipping file.")
        continue
    end
    if n_cells(real_spikes) < 5
        println("! Spikes for too few ($(n_cells(real_spikes))) cells. Skipping file.")
        continue
    end

    # create rasters

    # bin_size is now set at the commandline. Default value is 10ms
    println("  Computing spike rasters at bin size $(1000bin_size) ms")
    real_raster = raster(real_spikes, bin_size)
    sim_raster = raster(sim_spikes, bin_size)

    # gonna try only multiples of 5
    size_range = intersect(1:n_cells(real_spikes), default_size_range)
    println("  Preparing to fit models to subsamples of sizes $(join(size_range, ", "))")
    # From what I've seen in the JLD docs, there's no "append" mode. So first,
    # we have to load the existing data from the jld file, if it exists.

    distros = Dict{String,AbstractBinaryVectorDistribution}()
    if isfile(joinpath(information_dir, "$root_name-$rec_idx.jld"))
        print("  Found $root_name-$rec_idx.jld. Loading existing distributions.")
        try
            distros = load(joinpath(information_dir, "$root_name-$rec_idx.jld"))
        catch y
            println("! Exception occurred. Ignoring $root_name-$rec_idx.jld")
            distros = Dict{String,AbstractBinaryVectorDistribution}()
            show(y)
        end
    end
    try
        jldopen(joinpath(information_dir, "$root_name-$rec_idx.jld"), "w") do file
            for sample_size in size_range
                index_set = zeros(Int, sample_size)
                for trial = 1:min(n_trials, binomial(n_cells(real_spikes), sample_size))
                    sort!(random_subset!(1:n_cells(real_spikes), index_set))
                    III = index_set_to_int(index_set)
                    println("    size = $sample_size, trial $trial: [$(join(index_set,","))]")
                    print("      Writing to $root_name-$rec_idx.jld: ")
                    for XXX in ["1","2","N"]
                        for YYY in ["real","sim"]
                            distro_name = "P_$(XXX)_$(YYY)_$III"
                            if !haskey(distros, distro_name) || metadata(distros[distro_name], :CRCNS_script_version, v"0.1") < CRCNS_script_version
                                # sometimes, even though a conditional statement
                                # *must* define a variable, julia complains at me
                                # because it thinks that variable is not defined.
                                # So, just in case, I'm doing this dumb shit.
                                P = 0
                                if XXX == "1"
                                    P = first_order_model(YYY == "real" ? real_raster : sim_raster, index_set;
                                        CRCNS_script_version=CRCNS_script_version, verbose=verbose,
                                        source="CRCNS/$root_name-$rec_idx ($YYY)", bin_size=bin_size)
                                elseif XXX == "2"
                                    P = second_order_model(YYY == "real" ? real_raster : sim_raster, index_set;
                                        CRCNS_script_version=CRCNS_script_version, verbose=verbose,
                                        source="CRCNS/$root_name-$rec_idx ($YYY)", bin_size=bin_size, maxtime=maxtime)
                                else
                                    P = data_model(YYY == "real" ? real_raster : sim_raster, index_set;
                                        CRCNS_script_version=CRCNS_script_version, verbose=verbose,
                                        source="CRCNS/$root_name-$rec_idx ($YYY)", bin_size=bin_size)
                                end
                                if n_bits(P) <= Probability.ISING_METHOD_THRESHOLD
                                    entropy(P)
                                end
                                distros[distro_name] = P
                                write(file, distro_name, P)
                                print("$YYY/$XXX, ")
                            end
                        end
                    end
                    println()
                end
            end
        end
    catch y
        println()
        println("! Exception occurred while processing file $(joinpath(information_dir, "$root_name-$rec_idx.jld")). Skipping file.")
        show(y)
        continue
    end
    push!(successful_files, sim_file)
end

println("Successfully processed these files:")
println("\t$(join(successful_files,"\n\t"))")

println("CRCNS_output_information: END SCRIPT $(now())")
println("-" ^ 80)
