
using MAT

"""
    CRCNS_Stimulus(file_name, rec_idx; kwargs...)

A wrapper class for manipulating/computing with the stimulus used in the CRCNS
experiments.
"""
type CRCNS_Stimulus <: AbstractStimulus
    file_loc::String
    file_name::String
    rec_idx::Int # internally-used index for accessing array/dictionary elements
    rec_no::Int  # The actual recording number per the experimenter's labeling.
    bits::BitMatrix
    N::Tuple{Int,Int}
    px::Tuple{Int,Int}
    d::Tuple{Int,Int}
    mm_per_px::Float64
    frame_length_s::Float64
    onset::Float64
    metadata::Dict{Any,Any}

    function CRCNS_Stimulus(mat_file::String, recording_index::Int; verbose=false, kwargs...)
        fname = basename(mat_file)
        floc = dirname(abspath(mat_file))
        if verbose
            println("Reading CRCNS stimulus data from file: $(joinpath(floc, fname))")
        end
        vars = matread(joinpath(floc, fname))

        px = (Int(vars["stimulus"]["param"][recording_index]["x"]), Int(vars["stimulus"]["param"][recording_index]["y"]))
        d = (Int(vars["stimulus"]["param"][recording_index]["dx"]), Int(vars["stimulus"]["param"][recording_index]["dy"]))
        N = map(Int, px ./ d)
        mm_per_px = vars["stimulus"]["pixelsize"][recording_index] / 1000.0
        frame_length_s = vars["stimulus"]["frame"][recording_index]
        onset = vars["stimulus"]["onset"][recording_index]

        if verbose
            println("Loading binary data from $(joinpath(floc,"..","ran1.bin"))")
        end
        bin_file = joinpath(floc,"..","ran1.bin")

        # now load the frames. remember that read! will throw an error if the
        # number of bits to read is not 0 mod 64 and the bits beyond that are
        # nonzero. It's a bug with Julia but it's not too bad to work around.
        N_frames = Int(vars["stimulus"]["Nframes"][recording_index])
        temp_N_frames = N_frames
        while mod(prod(N) * temp_N_frames, 64) != 0
            temp_n_frames += 1
        end
        stimulus = falses(prod(N), temp_N_frames)
        read!(bin_file, stimulus)
        stimulus = stimulus[:,1:N_frames]

        new(floc, fname, recording_index, vars["datainfo"]["RecNo"][recording_index],
            stimulus, N, px, d, mm_per_px, frame_length_s, onset, Dict(kwargs))
    end
end

frame_size(S::CRCNS_Stimulus) = S.px
frame_time(S::CRCNS_Stimulus) = S.frame_length_s

matrix_form(S::CRCNS_Stimulus) = (-1.0) .^ (!S.bits)

function show(io::IO, S::CRCNS_Stimulus)
    println(io, "CRCNS binary whitenoise stimulus")
    println(io, "Used for file $(S.file_name), recording $(S.rec_no)")
    println(io, "Frame size (w,h): $(frame_size(S)) pixels, $(S.mm_per_px .* frame_size(S)) mm")
    println(io, "Framerate: $(frame_rate(S)) Hz ($(frame_time(S)) s/frame)")
end

"""
    frame_image(S::CRCNS_Stimulus, idx::Int)
    frame_image(S::CRCNS_Stimulus, t::Float64, relative_time=false)

Returns a matrix of size `reverse(S.px)` which is the actual image displayed on
screen. `idx` is the index of the frame, while `t` is the time of display. If
`relative_time` is `true`, time `t = 0.0` corresponds to stimulus onset.

"""
# TODO all of these will need checking to make sure I am consistent in reshaping stuff.
function frame_image(S::CRCNS_Stimulus, idx::Int)
    # first get the column from bits
    # frame_vec = S.bits[:,idx]
    frame = (-1.0) .^ (!reshape_frame_vec(S, S.bits[:,idx]))
    return kron(frame, ones(S.d))
end
frame_image(S::CRCNS_Stimulus, t::Float64, relative_time::Bool=false) = frame_image(S, ceil(Int, (relative_time ? t - S.onset : t)/frame_time(S)))
reshape_frame_vec(S::CRCNS_Stimulus, fv) = reshape(fv, reverse(S.N))

compute_STRFs(spike_hist::Matrix{Float64}, S::CRCNS_Stimulus; kwargs...) = _compute_STRFs(spike_hist, S; kwargs...)
