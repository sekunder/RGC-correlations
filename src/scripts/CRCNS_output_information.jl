
include("../probability/probability.jl")
include("../spikes/spikes.jl")
include("../stimulus/stimulus.jl")
using Spikes, Stimulus, Probability
include("../util/constants.jl")
include("../util/misc.jl")
include("../util/CRCNS/analysis_functions.jl")
using JLD

# plan: loop through folder of simulated STRFs (those files also have spike
# trains!):

# 1. load the simulated and real spike trains. compute rasters.

# 2. loop through: sample sizes, number of trials. Fit P_*, compute entropy (if
# appropriate), add this to a dictionary structure, etc.

sim_jld_dir = joinpath(CRCNS_STRF_dir, "sim")

sim_jld_files = filter(x -> endswith(x,".jld"), readdir(sim_jld_dir))

println("-" ^ 80)
println("CRCNS_output_information: BEGIN SCRIPT $(now())")

n_trials = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 20

for dir in [CRCNS_information_dir]
    if !isdir(dir)
        println("* Creating path $dir")
        mkpath(dir)
    end
end

println("Loading simulated and real spikes, then fitting P_1, P_2 to various subsamples of the data.")
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
            throw(y)
        end
    end
    if rec_idx == 0
        println("! Found $sim_file, but could not determine recording index.")
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
        println("! CRCNS spikes for $(n_cells(real_spikes)) cells, but simulated data for $(n_cells(sim_spikes)). Skipping...")
        continue
    end

    # create rasters
    #TODO maybe try a range of time bins? Most papers used 10ms, and frame_time tends to be ~16ms.
    sim_STRFs = load(joinpath(sim_jld_dir, sim_file), "STRFs")
    dt = frame_time(sim_STRFs[1]); sim_STRFs = 0 # hint for garbage collection
    real_raster = raster(real_spikes, dt)
    sim_raster = raster(sim_spikes, dt)

    # decide sample size ranges. I think 5:5:40 is about right, intersected with 1:N_neurons, then union N_neurons
    size_range = sort!(union(intersect(1:n_cells(real_spikes), 5:5:40), n_cells(real_spikes)))
    println("  Preparing to fit models to subsamples for sizes $(join(size_range, ", "))")
    # From what I've seen in the JLD docs, there's no "append" mode. So first,
    # we have to load the existing data from the jld file, if it exists.
    P_1_real = Dict{Int, BernoulliCodeDistribution}()
    P_2_real = Dict{Int, IsingDistribution}()
    P_N_real = Dict{Int, DataDistribution}()
    P_1_sim = Dict{Int, BernoulliCodeDistribution}()
    P_2_sim = Dict{Int, IsingDistribution}()
    P_N_sim = Dict{Int, DataDistribution}()
    if isfile(joinpath(CRCNS_information_dir, "$root_name.jld"))
        P_1_real = load(joinpath(CRCNS_information_dir, "$root_name.jld"), "P_1_real")
        P_2_real = load(joinpath(CRCNS_information_dir, "$root_name.jld"), "P_2_real")
        P_N_real = load(joinpath(CRCNS_information_dir, "$root_name.jld"), "P_N_real")
        P_1_sim = load(joinpath(CRCNS_information_dir, "$root_name.jld"), "P_1_sim")
        P_2_sim = load(joinpath(CRCNS_information_dir, "$root_name.jld"), "P_2_sim")
        P_N_sim = load(joinpath(CRCNS_information_dir, "$root_name.jld"), "P_N_sim")
    end
    jldopen(joinpath(CRCNS_information_dir, "$root_name.jld"), "w") do file
        #TODO loop through size_range, loop for n_trials, etc.
        for sample_size in size_range
            #TODO add some verbosity
            
            # TODO check what subsets have already been fit with at least the
            # current version of the script

            # TODO but in order to do that I'll have to include a version in the
            # metadata of all the objects? Well, not necessarily.
            index_set = zeros(Int, sample_size)
            for trial = 1:min(n_trials, binomial(n_cells(real_spikes), sample_size))
                random_subset!(1:n_cells(real_spikes), index_set)
                index_int = index_set_to_int(index_set)
                if !haskey(P_1_real, index_int) || metadata(P_1_real[index_int], :CRCNS_script_version, v"0.1") < CRCNS_script_version
                    # short-circuit operator!
                    P_1_real[index_int] = first_order_model(real_raster, index_set; CRCNS_script_version=CRCNS_script_version)
                end
                if !haskey(P_1_sim, index_int) || metadata(P_1_sim[index_int], :CRCNS_script_version, v"0.1") < CRCNS_script_version
                    P_1_sim[index_int] = first_order_model(sim_raster, index_set; CRCNS_script_version=CRCNS_script_version)
                end
                if !haskey(P_2_real, index_int) || metadata(P_2_real[index_int], :CRCNS_script_version, v"0.1") < CRCNS_script_version
                    # short-circuit operator!
                    P_2_real[index_int] = second_order_model(real_raster, index_set; CRCNS_script_version=CRCNS_script_version)
                end
                if !haskey(P_2_sim, index_int) || metadata(P_2_sim[index_int], :CRCNS_script_version, v"0.1") < CRCNS_script_version
                    P_2_sim[index_int] = second_order_model(sim_raster, index_set; CRCNS_script_version=CRCNS_script_version)
                end
                if !haskey(P_N_real, index_int) || metadata(P_N_real[index_int], :CRCNS_script_version, v"0.1") < CRCNS_script_version
                    # short-circuit operator!
                    P_N_real[index_int] = data_model(real_raster, index_set; CRCNS_script_version=CRCNS_script_version)
                end
                if !haskey(P_N_sim, index_int) || metadata(P_N_sim[index_int], :CRCNS_script_version, v"0.1") < CRCNS_script_version
                    P_N_sim[index_int] = data_model(sim_raster, index_set; CRCNS_script_version=CRCNS_script_version)
                end
            end
        end
        file["P_1_real"] = P_1_real
        file["P_2_real"] = P_2_real
        file["P_N_real"] = P_N_real
        file["P_1_sim"] = P_1_sim
        file["P_2_sim"] = P_2_sim
        file["P_N_sim"] = P_N_sim
    end
end

println("CRCNS_output_information: END SCRIPT $(now())")
println("-" ^ 80)
