#TODO(medium) make a CRCNS folder in /src/stimulus
using MAT
"""
    CRCNS_Stimulus(mat_file, recording_index;
        verbose=false, kwargs...)

Returns a `GrayScaleStimulus` object which represents the stimulus used in the
indicated recording from the given file. Includes this information in the
metadata.

"""
function CRCNS_Stimulus(mat_file::String, recording_index::Int; verbose=false, kwargs...)
    fname = basename(mat_file)
    floc = dirname(abspath(mat_file))
    if verbose
        #TODO(medium) Fix all verbose output to include function name.
        println("Reading CRCNS stimulus data from file: $(joinpath(floc, fname))")
    end
    vars = matread(joinpath(floc, fname))

    px = (Int(vars["stimulus"]["param"][recording_index]["x"]), Int(vars["stimulus"]["param"][recording_index]["y"]))
    d = (Int(vars["stimulus"]["param"][recording_index]["dx"]), Int(vars["stimulus"]["param"][recording_index]["dy"]))
    N = map(Int, px ./ d)
    mm_per_px = vars["stimulus"]["pixelsize"][recording_index] / 1000.0
    frame_length_s = vars["stimulus"]["frame"][recording_index]
    frame_rate = 1.0 / frame_length_s
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
        temp_N_frames += 1
    end
    stimulus = falses(prod(N), temp_N_frames)
    read!(bin_file, stimulus)
    stimulus = stimulus[:,1:N_frames]

    # new(floc, fname, recording_index, vars["datainfo"]["RecNo"][recording_index],
    #     stimulus, N, px, d, mm_per_px, frame_length_s, onset, Dict(kwargs))

    # create some appropriate metadata
    meta = Dict()
    meta[:file] = joinpath(floc,fname)
    meta[:recording] = "Recording $(vars["datainfo"]["RecNo"][recording_index]) (recording index $(recording_index))"
    meta[:source] = "CRCNS Binary Whitenoise stimulus"
    return GrayScaleStimulus(stimulus, N, px, d, mm_per_px, frame_length_s, frame_rate, onset, true, merge(Dict(kwargs),meta))
end
