
using MAT
"""
    CRCNS_Stimulus(mat_file, recording_index;
        verbose=0, kwargs...)

Returns a `GrayScaleStimulus` object which represents the stimulus used in the
indicated recording from the given file. Includes this information in the
metadata.

"""
@everywhere function CRCNS_Stimulus(mat_file::String, recording_index::Int; verbose=0, single_rec=false, kwargs...)
    fname = basename(mat_file)
    floc = dirname(abspath(mat_file))
    if verbose > 0
        println("CRCNS_Stimulus: Reading CRCNS stimulus data from file: $(joinpath(floc, fname))")
    end
    vars = matread(joinpath(floc, fname))

    # because matlab is terrible (or perhaps just the way MAT opens matlab files
    # is terrible...), I need to work around the possibility that the number of
    # recording indexes is 1, meaning that the intermediate [recording_index] in
    # the expressions below will actually throw an error.
    if single_rec
        px = [Int(vars["stimulus"]["param"]["x"]), Int(vars["stimulus"]["param"]["y"])]
        d = [Int(vars["stimulus"]["param"]["dx"]), Int(vars["stimulus"]["param"]["dy"])]
        N = map(Int, px ./ d)
        mm_per_px = vars["stimulus"]["pixelsize"] / 1000.0
        frame_length_s = vars["stimulus"]["frame"]
        frame_rate = 1.0 / frame_length_s
        onset = vars["stimulus"]["onset"]
    else
        px = [Int(vars["stimulus"]["param"][recording_index]["x"]), Int(vars["stimulus"]["param"][recording_index]["y"])]
        d = [Int(vars["stimulus"]["param"][recording_index]["dx"]), Int(vars["stimulus"]["param"][recording_index]["dy"])]
        N = map(Int, px ./ d)
        mm_per_px = vars["stimulus"]["pixelsize"][recording_index] / 1000.0
        frame_length_s = vars["stimulus"]["frame"][recording_index]
        frame_rate = 1.0 / frame_length_s
        onset = vars["stimulus"]["onset"][recording_index]
    end

    if verbose > 0
        println("CRCNS_Stimulus: Loading binary data from $(joinpath(floc,"..","ran1.bin"))")
    end
    bin_file = joinpath(floc,"..","ran1.bin")

    # now load the frames. remember that read! will throw an error if the
    # number of bits to read is not 0 mod 64 and the bits beyond that are
    # nonzero. It's a bug with Julia but it's not too bad to work around.
    if single_rec
        N_frames = Int(vars["stimulus"]["Nframes"])
    else
        N_frames = Int(vars["stimulus"]["Nframes"][recording_index])
    end
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
    RecNo = Int(single_rec ? vars["datainfo"]["RecNo"] : vars["datainfo"]["RecNo"][recording_index])
    meta[:recording] = "Recording $RecNo (recording index $(recording_index))"
    meta[:source] = "CRCNS Binary Whitenoise stimulus"
    meta[:mm_per_px] = mm_per_px
    return GrayScaleStimulus(stimulus, N, px, d, frame_length_s, onset, true, merge(Dict(kwargs),meta))
end
