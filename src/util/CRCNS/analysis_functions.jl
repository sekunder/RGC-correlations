# miscellaneous functions that don't quite fit in any module, since they might
# use features from multiple modules.
using JLD
"""
    CRCNS_output_STRFs(mat_file, output_dir, rec_idx; kwargs...)

Comutes the spatio-temporal receptive fields from the data available in
`mat_file`, returning the result as an array of `GrayScaleStimulus` objects.
Moreover, saves this array to `output_dir` using the `JLD` package, together
with some extra information (a time stamp and a version number, so I know which
version of the `read`-type functions to use.)

"""
function CRCNS_output_STRFs(mat_file, rec_idx, output_dir=dirname(abspath(mat_file));
    verbose=false, CRCNS_script_version=v"0.1", kwargs...)
    # fname = basename(mat_file)
    # floc = dirname(abspath(mat_file))
    if !isdir(output_dir)
        if verbose
            println("CRCNS_output_STRFs: making directory $(abspath(output_dir))")
        end
        mkpath(output_dir)
    end

    stim = CRCNS_Stimulus(mat_file, rec_idx; verbose=verbose)
    spikes = Spikes.CRCNS_get_spikes_from_file(mat_file, rec_idx)
    idx = mapreduce(x -> 2^(x-1), +, spikes.I)

    spike_hist = histogram(spikes, frame_time(stim))
    STRFs = compute_STRFs(spike_hist, stim)
    timestamp = now()
    if verbose
        println("CRCNS_output_STRFs: writing jld file to $output_dir")
    end
    save(joinpath(output_dir, "CRCNS_STRF_$idx.jld"), "CRCNS_script_version", CRCNS_script_version, "STRFs", STRFs, "timestamp", timestamp)
    return STRFs
end
