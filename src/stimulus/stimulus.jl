"""
Stimulus

A wrapper class and API for handling stimuli without needing to worry about the
specific implementation for that experiment.
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

function compute_STRFs(spike_raster::Union{Matrix{Bool}, BitMatrix}, stim::AbstractStimulus; kwargs...)
    # TODO implement this. The idea is that frame_vector will return a matrix
    # where each column is one frame of the stimulus, so this function should
    # generically work on any stimulus with frame_vector implemented.
end

################################################################################
#### Include definitions of concrete types
################################################################################
include("CRCNS_Stimulus.jl")

export frame_matrix, frame_vector, frame_image,
       CRCNS_Stimulus
end
