
using MAT
"""
    Ganmor_Stimulus(file; verbose=0, kwargs...)

Returns a `GrayScaleStimulus` object for the video used in the Ganmor
experiments.

"""
function Ganmor_Stimulus(filename; verbose=0, kwargs...)

    movie = matread(filename)["Movie"]
    # movie = read(filename, "Movie") # MAT.jl has some weird bug, not gonna try and figure it out.
    N_frames, ht, wd = size(movie)
    N = [wd, ht]
    d = [1,1]
    px = N .* d

    mm_per_px = 1.0
    frame_rate = 30#TODO look at the paper
    frame_length_s = 1.0 /frame_rate
    zerotonegative = true
    onset = 0.0

    px_vals = zeros(UInt8, ht * wd, N_frames)
    # TODO problem: the video, if converted to float64, will be ~9gb in size.
    # ostensibly, this is not the end of the world, since i'll need to perform
    # two expensive calculations with it, after which I'll have much more
    # manageable objects to deal with. but still, how should I handle this?
    for frame = 1:N_frames
        px_vals[:,frame] = movie[frame,:,:][:]
    end

    meta = Dict()
    meta[:file] = filename
    meta[:source] = "Ganmor data, Natural scene movie"

    return GrayScaleStimulus(px_vals, N, px, d, mm_per_px, frame_length_s, frame_rate, onset, zerotonegative, merge(Dict(kwargs), meta))
end
