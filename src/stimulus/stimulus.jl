"""
Stimulus

A wrapper class and API for handling stimuli without needing to worry about the
specific implementation for that experiment. Basically, I want to be able to
call `compute_STRFs(raster, stimulus)` regardless of the data source.
"""
module Stimulus

import Base.show

################################################################################
#### Abstract type to use as a placeholder/to get all docstrings
################################################################################
""""
    AbstractStimulus

Abstract type used as placeholder.
"""
abstract AbstractStimulus

################################################################################
#### Miscellaneous functions and constants
################################################################################
include("../util/metadata.jl")

default_CRCNS_dir = ispath("/Users/sekunder/Documents/MATLAB/crcns_ret-1") ? "/Users/sekunder/Documents/MATLAB/crcns_ret-1" : "/data1/homes/abk170/crcns_ret-1"

"""
    framerate(S::AbstractStimulus)

Convenience method returns `1/frametime(S)`. Each concrete subtype of
AbstractStimulus should implement `frametime`.
"""
frame_rate(S::AbstractStimulus) = 1.0 / frame_time(S)

"""
    compute_STRFs(spike_hist, stim; kwargs...)

Generic implementation of STRF computation. Uses `matrix_form(stim)` to perform
computations. Returns a matrix.

Life will be easy if we assume that the bins of the histogram are 1-to-1 with
the frames of the stimulus. So, implement concrete subtypes of
`AbstractStimulus` with that in mind.

For keyword arguments: specifying `window_length_frames` overrides any value
computed from `window_length_s`
"""
function _compute_STRFs(spike_hist::Matrix{Float64}, stim::AbstractStimulus; kwargs...)
    # preprocessing stuff: get the window length
    dkwargs = Dict(kwargs)
    window_length_s = get(dkwargs, :window_length_s, 0.5)
    window_length_frames = round(Int, window_length_s / frame_time(stim))
    if haskey(dkwargs, :window_length_frames)
        window_length_frames = dkwargs[:window_length_frames]
    end

    # get the matrix form of the stimulus
    stimulus = matrix_form(stim)

    N_bins, N_cells = size(spike_hist)
    N_bits_per_frame, N_frames = size(stimulus)
    if N_bins != N_frames
        error("_compute_STRFs: mismatched histogram and matrix form (N_bins = $N_bins, N_frames = $N_frames)")
    end

    avg_spikes_per_frame = zeros(spike_hist)
    for j = 1:N_cells
        avg_spikes_per_frame[window_length_frames:end,j] = spike_hist[window_length_frames:end,j] / sum(spike_hist[window_length_frames:end-1,j])
    end

    RFs = zeros(N_bits_per_frame, N_cells, window_length_frames)
    for k = 1:window_length_frames
        RFs[:, :, window_length_frames - k + 1] = stimulus * avg_spikes_per_frame
        avg_spikes_per_frame = circshift(avg_spikes_per_frame, [-1,0])
    end
    return RFs
end

################################################################################
#### Include definitions of concrete types
################################################################################
include("GrayScaleStimulus.jl")
include("CRCNS_Stimulus.jl")

export show, frame_size, frame_time, frame_rate, matrix_form, frame_image,
       GrayScaleStimulus, CRCNS_Stimulus,
       compute_STRFs, STRF_response

end
