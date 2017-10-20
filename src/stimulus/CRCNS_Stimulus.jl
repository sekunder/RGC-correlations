"""
    CRCNS_Stimulus(filename, rec_idx; kwargs...)

A wrapper class for manipulating/computing with the stimulus used in the CRCNS
experiments.
"""
type CRCNS_Stimulus <: AbstractStimulus
    Bits::BitMatrix
    # TODO more fields go here


    function CRCNS_Stimulus(filename::String, recording_index::Int; kwargs...)
        # TODO write this function, and decide what fields this type needs
    end
end

function show(io::IO, S::CRCNS_Stimulus)
    # TODO display useful information.
end

function frame_matrix(S::CRCNS_Stimulus, I)
    # TODO return a 3d array (3rd dim might be 1, or omitted?)
end

function frame_vector(S::CRCNS_Stimulus, I)
    # TODO return a 2d array (might be a vector)
end

function frame_image(S::CRCNS_Stimulus, I)
    # TODO return a 3d array (3rd dim might be 1?)
end
