

include("../util/init.jl")
include("../util/CRCNS/analysis_functions.jl")
using JLD

# plan: loop through folder of simulated STRFs (those files also have spike
# trains!):

# 1. load the simulated and real spike trains. compute rasters.

# 2. loop through: sample sizes, number of trials. Fit P_*, compute entropy (if
# appropriate), add this to a dictionary structure, etc.

cline_args = process_args(ARGS, parse_defaults=Dict("n_trials"=>20, "bin_size"=>10e-3))
n_trials = cline_args["n_trials"]
# To make life easy, I want to be able to say julia /script 10 to do 10ms
dt = isa(cline_args["bin_size"], Integer) ? cline_args["bin_size"] / 1000 : cline_args["bin_size"]
verbose = cline_args["v"] ? 1 : cline_args["verbose"]

sim_jld_dir = joinpath(CRCNS_STRF_dir, "sim")

sim_jld_files = filter(x -> endswith(x,".jld"), readdir(sim_jld_dir))

println("-" ^ 80)
println("CRCNS_output_information: BEGIN SCRIPT $(now())")

for dir in [CRCNS_information_dir]
    if !isdir(dir)
        println("* Creating path $dir")
        mkpath(dir)
    end
end

println("Loading simulated and real spikes, then fitting P_1, P_2 to subsamples of the data.")
println("Will run $n_trials trials for each sample size")
println("Will write output to $CRCNS_information_dir")


for sim_file in sim_jld_files
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
    if !isfile(joinpath(CRCNS_data_dir, "$root_name.mat"))
        println("! Found $sim_file, but could not find $root_name.mat. Skipping file.")
        continue
    end
    println("* Processing $root_name, recording index $rec_idx")
    println("  Loading CRCNS spikes from $root_name.mat")
    real_spikes = Spikes.CRCNS_get_spikes_from_file(joinpath(CRCNS_data_dir, "$root_name.mat"), rec_idx)
    println("  Loading simulated spikes from $sim_file")
    sim_spikes = load(joinpath(sim_jld_dir, sim_file), "spikes")
    if n_cells(real_spikes) != n_cells(sim_spikes)
        println("! CRCNS spikes for $(n_cells(real_spikes)) cells, but simulated data for $(n_cells(sim_spikes)). Skipping file.")
        continue
    end
    if n_cells(real_spikes) < 5
        println("! Spikes for too few ($(n_cells(real_spikes))) cells. Skipping file.")
        continue
    end

    # create rasters

    # dt is now set at the commandline. Default value is 10ms
    println("  Computing spike rasters at bin size $(1000dt) ms")
    real_raster = raster(real_spikes, dt)
    sim_raster = raster(sim_spikes, dt)

    # gonna try only multiples of 5
    size_range = intersect(1:n_cells(real_spikes), 5:5:40)
    println("  Preparing to fit models to subsamples of sizes $(join(size_range, ", "))")
    # From what I've seen in the JLD docs, there's no "append" mode. So first,
    # we have to load the existing data from the jld file, if it exists.

    distros = Dict{String,AbstractBinaryVectorDistribution}()
    if isfile(joinpath(CRCNS_information_dir, "$root_name-$rec_idx.jld"))
        print("  Found $root_name-$rec_idx.jld. Loading existing distributions...")
        distros = load(joinpath(CRCNS_information_dir, "$root_name-$rec_idx.jld"))
        println("done.")
    end
    # distro_names = ["P_1_real", "P_2_real", "P_N_real", "P_1_sim", "P_2_sim", "P_N_sim"]
    # distros = [P_1_real, P_2_real, P_N_real, P_1_sim, P_2_sim, P_N_sim]
    jldopen(joinpath(CRCNS_information_dir, "$root_name-$rec_idx.jld"), "w") do file
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
                                    source="CRCNS/$root_name-$rec_idx ($YYY)", bin_size=dt)
                            elseif XXX == "2"
                                P = second_order_model(YYY == "real" ? real_raster : sim_raster, index_set;
                                    CRCNS_script_version=CRCNS_script_version, verbose=verbose,
                                    source="CRCNS/$root_name-$rec_idx ($YYY)", bin_size=dt)
                            else
                                P = data_model(YYY == "real" ? real_raster : sim_raster, index_set;
                                    CRCNS_script_version=CRCNS_script_version, verbose=verbose,
                                    source="CRCNS/$root_name-$rec_idx ($YYY)", bin_size=dt)
                            end
                            distros[distro_name] = P
                            write(file, distro_name, P)
                            print("$YYY/$XXX,")
                        end
                    end
                end
                println()
            end
        end
        # println("  Writing to file $(joinpath(CRCNS_information_dir, "$root_name-$rec_idx.jld"))")
        # for (n, v) in zip(distro_names, distros)
        #     write(file, n, v)
        # end
    end
end

println("CRCNS_output_information: END SCRIPT $(now())")
println("-" ^ 80)
