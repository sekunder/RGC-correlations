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
    verbose=0, CRCNS_script_version=v"0.1", fp=STDOUT, indent=0, kwargs...)
    # fname = basename(mat_file)
    # floc = dirname(abspath(mat_file))
    dkwargs = Dict(kwargs)

    base_name = remove_extension(basename(mat_file))
    if !isdir(output_dir)
        if verbose > 0
            println(fp, "$(ts() * sp(indent)) CRCNS_output_STRFs: making directory $(abspath(output_dir))")
        end
        mkpath(output_dir)
    end

    single_rec = pop!(dkwargs,:single_rec,false)
    # load spikes and stimulus from files
    stim = CRCNS_Stimulus(mat_file, rec_idx; verbose=verbose, single_rec=single_rec) # don't save stim, it's like ~30MB
    spikes = CRCNS_get_spikes_from_file(mat_file, rec_idx)

    # compute histogram and STRFs
    spike_hist = histogram(spikes, frame_time(stim); N_bins=n_frames(stim))
    STRFs = compute_STRFs(spike_hist, stim)

    # add all the useful metadata
    spikes.metadata[:frame_hist] = spike_hist
    hide_metadata!(spikes, :frame_hist)
    for rf in STRFs
        rf.metadata[:animal] = spikes.metadata[:animal]
    end

    # Trying a new approach to file io
    savespikes(spikes)
    map(savestimulus, STRFs)
    # idx = index_set_to_int(spikes.I)
    # filename = "$base_name-$(rec_idx)_STRF_$idx.jld"
    # #MAYBEDO figure out how to use expressions to make this customizable.
    #
    # spike_hist = histogram(spikes, frame_time(stim); N_bins=n_frames(stim))
    # # now let's check if these things have already been computed, to save time
    # file_exists = ispath(joinpath(output_dir, filename))
    # STRF_exists = false
    # if file_exists
    #     STRF_exists = try
    #             d = load(joinpath(output_dir, filename))
    #             d["CRCNS_script_version"] == CRCNS_script_version
    #         catch y
    #             if verbose > 0
    #                 println(fp, "$(ts() * sp(indent)) CRCNS_output_STRFs: error encountered while reading ")
    #             end
    #             false
    #         end
    # end
    # if !file_exists || !STRF_exists
    #     STRFs = compute_STRFs(spike_hist, stim)
    #     timestamp = now()
    # else
    #     if verbose > 0
    #         println(fp, "$(ts() * sp(indent)) CRCNS_output_STRFs: reading jld file $(joinpath(output_dir, filename))")
    #     end
    #     STRFs = d["STRFs"]
    # end
    # if !file_exists
    #     if verbose > 0
    #         println(fp, "$(ts() * sp(indent)) CRCNS_output_STRFs: writing jld file $(joinpath(output_dir, filename))")
    #     end
    #     save(joinpath(output_dir, filename), "CRCNS_script_version", CRCNS_script_version, "STRFs", STRFs, "timestamp", timestamp)
    # end
    return stim, spikes, spike_hist, STRFs
end

"""
    CRCNS_collect_I_2(dir=CRCNS_information_dir; kwargs...)

Convenience function. Loop through the given directory, opening all `.jld` files
it encounters that contain objects from the `Probability` module. For each index
set `I` that has a value of `P_*_real_I` and `P_*_sim_I` (for `* = 1, 2, N`),
computes the information ratio of the distributions, defined as
```latex
I_2 = (H[P_1] - H[P_2]) / (H[P_1] - H[P_N])
```

Returns `(I_2_real, I_2_sim)` where each is a `Dict{Int,Vector{Float64}}` with
`I_2_*[k]` the values of `I_2_*` for subsets of size `k`. Stores the values in
the same order (so `I_2_real[k][j]` is computed from the same index set as
`I_2_sim[k][j]`)

"""
function CRCNS_collect_entropy(dir=CRCNS_information_dir;
    verbose=0, CRCNS_script_version=CRCNS_script_version, fp=STDOUT, indent=0, kwargs...)

    jld_files = filter(x -> endswith(x, ".jld"), readdir(dir))
    if verbose > 0
        println(fp, "$(ts() * sp(indent)) CRCNS_collect_entropy: found $(length(jld_files)) files in $dir")
    end
    # H_all = Dict{Int,Matrix{Float64}}()
    I_2_real = Dict{Int,Vector{Float64}}()
    I_2_sim = Dict{Int,Vector{Float64}}()

    for filename in jld_files
        if verbose > 0
            println(fp, "$(ts() * sp(indent)) CRCNS_collect_entropy: Inspecting file $filename")
        end
        H_temp = Dict()
        try
            distros = load(joinpath(dir,filename))
            for (dname, P) in distros
                PPP,XXX,YYY,III = split(dname,"_")
                if !haskey(H_temp,XXX)
                    H_temp[XXX] = Dict()
                end
                if !haskey(H_temp[XXX], YYY)
                    H_temp[XXX][YYY] = Dict()
                end
                H_temp[XXX][YYY][III] = entropy(P)
            end
        catch y
            if verbose > 0
                println(fp, "$(ts() * sp(indent+1)) Exception occurred. Skipping file.")
                show(y)
                continue
            end
        end
        valid_Is = intersect([keys(H[XXX][YYY]) for XXX in ["1","2","N"], YYY in ["real","sim"]]...)
        for III in valid_Is
            k = count_ones(III)
            # H_all[k] = hcat(get!(H_all,k,zeros(6,0)), [H_temp["1"]["real"][III], H_temp["2"]["real"][III], H_temp["N"]["real"][III], H_temp["1"]["sim"][III], H_temp["2"]["sim"][III], H_temp["N"]["sim"][III]])
            # H_all[k] = hcat(get!(H_all,k,zeros(6,0)), [H_temp[XXX][YYY][III] for XXX in ["1","2","N"], YYY in ["real","sim"]])
            i2r = (H["1"]["real"][III] - H["2"]["real"][III]) / (H["1"]["real"][III] - H["N"]["real"][III])
            i2s = (H["1"]["sim"][III] - H["2"]["sim"][III]) / (H["1"]["sim"][III] - H["N"]["sim"][III])
            if isreal([i2r, i2s]) && all(isfinite.([i2r, i2s]))
                push!(get!(I_2_real, k, Float64[]), i2r)
                push!(get!(I_2_sim, k, Float64[]), i2s)
            else
                if verbose > 0
                    println(fp, "$(ts() * sp(indent+1)) Found bad index $III (size $k):")
                    println(fp, "$(ts() * sp(indent+2)) I_2_real = $i2r")
                    println(fp, "$(ts() * sp(indent+2)) I_2_sim  = $i2s")
                end
            end
        end
    end
    return I_2_real, I_2_sim
end
#OLD DOCSTRING:  returns `(H_1_real, H_2_real,
# H_N_real, H_1_sim, H_2_sim, H_N_sim)` where each is a
# `Dict{Int,Vector{Float64}}`, `H_*[k]` is the list of entropy values for
# distributions on `k` neurons.

# """
#     estimate_RF(strf::GrayScaleStimulus)
#
# Given an STRF, computes the variance in the time dimension; filters out
# """
# function estimate_RF(strf::GrayScaleStimulus)
# end

"""
    CRCNS_strf_stats(rf; zscore=3)

Fit a 2d guassian to the "variance image" of the given STRF. That is, first reduce the
temporal dimension by variance, then set everything that is less than `zscore` deviations
away from the mean to 0. Fit a gaussian to the resulting full image (i.e., the image is
scaled up to `frame_size(rf)`). Returns the center and covariance matrix of the resulting
gaussian in "pixel space."

If none of the pixels vary enough, returns `missing, missing`.

"""
function CRCNS_strf_stats(rf::GrayScaleStimulus; zscore=3)
    M = matrix_form(rf)
    V = var(M, 2)
    Z = (V - mean(V)) / sqrt(var(V))
    V[abs.(Z) .< zscore] = 0.0
    if all(V .== 0)
        return missing, missing
    end
    img = kron(reshape(V, resolution(rf)...), ones(patch_size(rf)...))
    C = imagemean(img)
    Σ = matrixcovariance(img)
    return C,Σ
end


for r in eachrow(dfstrf_real)
    rf = loadstimulus(r[:hash])
    r[:mm_per_px] = rf.mm_per_px
    C,S = CRCNS_strf_stats(rf)
    if !ismissing(S)
        va,ve = eig(S)
        r[:center] = C
        r[:v1] = va[1] * ve[:,1]
        r[:v2] = va[2] * ve[:,2]
    end
end
