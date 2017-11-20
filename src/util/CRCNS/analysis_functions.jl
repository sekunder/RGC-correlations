# miscellaneous functions that don't quite fit in any module, since they might
# use features from multiple modules.

using JLD
# include("../misc.jl")
"""
    CRCNS_output_STRFs(mat_file, output_dir, rec_idx; kwargs...)

Convenience function; performs a few related operations at once. Given a
filename and index, loads the appropriate spike trains and stimulus, then
computes the STRFs. Returns the stimulus, spikes, spike histogram, and array of
STRFs computed.

Moreover, saves this array to `output_dir` using the `JLD` package, together
with some extra information (a time stamp and a version number, so I know which
version of the `read`-type functions to use.)

Keyword argument `verbose` can be set to `0`,`1`, or `2`; `0` means no output;
`1` means output just from this function (i.e. it prints where the files are
being saved), `2` means verbose output from functions called within this one.

"""
function CRCNS_output_STRFs(mat_file, rec_idx, output_dir=dirname(abspath(mat_file));
    verbose=0, CRCNS_script_version=v"0.1", kwargs...)
    # fname = basename(mat_file)
    # floc = dirname(abspath(mat_file))
    dkwargs = Dict(kwargs)

    base_name = remove_extension(basename(mat_file))
    if !isdir(output_dir)
        if verbose > 0
            println("CRCNS_output_STRFs: making directory $(abspath(output_dir))")
        end
        mkpath(output_dir)
    end

    single_rec = pop!(dkwargs,:single_rec,false)
    stim = CRCNS_Stimulus(mat_file, rec_idx; verbose=verbose, single_rec=single_rec)
    spikes = Spikes.CRCNS_get_spikes_from_file(mat_file, rec_idx)
    idx = index_set_to_int(spikes.I)
    filename = "$base_name-$(rec_idx)_STRF_$idx.jld"
    #MAYBEDO figure out how to use expressions to make this customizable.

    spike_hist = histogram(spikes, frame_time(stim); N_bins=n_frames(stim))
    # now let's check if these things have already been computed, to save time
    file_exists = ispath(joinpath(output_dir, filename))
    STRF_exists = false
    if file_exists
        d = load(joinpath(output_dir, filename))
        STRF_exists = d["CRCNS_script_version"] == CRCNS_script_version
    end
    if !file_exists || !STRF_exists
        STRFs = compute_STRFs(spike_hist, stim)
        timestamp = now()
    else
        if verbose > 0
            println("CRCNS_output_STRFs: reading jld file $(joinpath(output_dir, filename))")
        end
        STRFs = d["STRFs"]
    end
    if !file_exists
        if verbose > 0
            println("CRCNS_output_STRFs: writing jld file $(joinpath(output_dir, filename))")
        end
        save(joinpath(output_dir, filename), "CRCNS_script_version", CRCNS_script_version, "STRFs", STRFs, "timestamp", timestamp)
    end
    return stim, spikes, spike_hist, STRFs
end

"""
    CRCNS_collect_entropy(dir=CRCNS_information_dir; kwargs...)

Convenience function. Loop through the given directory, opening all `.jld` files
it encounters that contain objects from the `Probability` module. For each index
set that has a value of `P_*_real` and `P_*_sim` (for `* = 1, 2, N`), it
computes the entropy of each distribution, then returns `H_all`, a
`Dict{Int,Matrix{Float64}}` where `H_all[k]` has six rows and some number of
columns. The rows are, in order, the entropy of `P_1_real`, `P_2_real`,
`P_N_real`, `P_1_sim`, `P_2_sim`, `P_N_sim`. Each column represents one sample
from one file.

"""
function CRCNS_collect_entropy(dir=CRCNS_information_dir;
    verbose=0, CRCNS_script_version=v"0.2", kwargs...)

    jld_files = filter(x -> endswith(x, ".jld"), readdir(dir))
    if verbose > 0
        println("CRCNS_collect_entropy: found $(length(jld_files)) files in $dir")
    end
    # H_1_real = Dict{Int,Vector{Float64}}()
    # H_2_real = Dict{Int,Vector{Float64}}()
    # H_N_real = Dict{Int,Vector{Float64}}()
    # H_1_sim = Dict{Int,Vector{Float64}}()
    # H_2_sim = Dict{Int,Vector{Float64}}()
    # H_N_sim = Dict{Int,Vector{Float64}}()
    H_all = Dict{Int,Matrix{Float64}}()

    distro_names = ["P_1_real", "P_2_real", "P_N_real", "P_1_sim", "P_2_sim", "P_N_sim"]

    for filename in jld_files
        if verbose > 0
            println("CRCNS_collect_entropy: Inspecting file $filename")
        end
        try
            distros = load(joinpath(dir,filename))
            if haskeys(distros, distro_names...)
                if verbose > 0
                    print("\tFound distributions. Computing entropies...")
                end
                # reminder: the function I keep forgetting is count_ones.
                index_ints = intersect([keys(distros[dname]) for dname in distro_names]...)
                for index_int in index_ints
                    k = count_ones(index_int)

                    if k <= Probability.ISING_METHOD_THRESHOLD
                        # push!(get!(H_1_real, k, Float64[]), entropy(distros["P_1_real"][index_int]))
                        # push!(get!(H_2_real, k, Float64[]), entropy(distros["P_2_real"][index_int]))
                        # push!(get!(H_N_real, k, Float64[]), entropy(distros["P_N_real"][index_int]))
                        # push!(get!(H_1_sim, k, Float64[]), entropy(distros["P_1_sim"][index_int]))
                        # push!(get!(H_2_sim, k, Float64[]), entropy(distros["P_1_sim"][index_int]))
                        # push!(get!(H_N_sim, k, Float64[]), entropy(distros["P_N_sim"][index_int]))
                        # ents = [entropy(distros["P_1_real"][index_int]), entropy(distros["P_2_real"][index_int]), entropy(distros["P_N_real"][index_int]), entropy(distros["P_1_sim"][index_int]), entropy(distros["P_1_sim"][index_int]), entropy(distros["P_N_sim"][index_int])]
                        ents = map(entropy, [distros[dn][index_int] for dn in distro_names])
                        if !haskey(H_all, k)
                            H_all[k] = zeros(6,0)
                        end
                        H_all[k] = hcat(H_all[k], ents)
                        if verbose > 1
                            print("$index_int,")
                        end
                    elseif verbose > 1
                        print("($index_int),")
                    end
                end
                if verbose > 0
                    println("done.")
                end
            elseif verbose > 0
                # couldn't find all necessary distros
                println("\tMissing keys $(join(setdiff(distro_names, keys(distros)),", ")), skipping file")
            end
            # if we've made it this far, I think I have not changed any of the
            # dictionaries, just modified the probability distributions
            # (entropy() has the side effect of storing the value as metadata).
            # So, I should be able to rewrite the file
            if verbose > 0
                print("\tWriting distros (with entropies) to $filename...")
            end
            jldopen(joinpath(dir,filename),"w") do file
                for (dist_name, dist) in distros
                    write(file, dist_name, dist)
                    if verbose > 1
                        print("$dist_name ")
                    end
                end
            end
            if verbose > 0
                println("done")
            end
        catch y
            if verbose > 0
                println()
                println("\tError: $y")
                println("\tskipping file")
            end
        end
    end
    # return H_1_real, H_2_real, H_N_real, H_1_sim, H_2_sim, H_N_sim
    return H_all
end
#OLD DOCSTRING:  returns `(H_1_real, H_2_real,
# H_N_real, H_1_sim, H_2_sim, H_N_sim)` where each is a
# `Dict{Int,Vector{Float64}}`, `H_*[k]` is the list of entropy values for
# distributions on `k` neurons.
