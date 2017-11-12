# miscellaneous functions that don't quite fit in any module, since they might
# use features from multiple modules.

using JLD
include("../misc.jl")
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
    base_name = remove_extension(basename(mat_file))
    if !isdir(output_dir)
        if verbose > 0
            println("CRCNS_output_STRFs: making directory $(abspath(output_dir))")
        end
        mkpath(output_dir)
    end

    stim = CRCNS_Stimulus(mat_file, rec_idx; verbose=(verbose > 1))
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
