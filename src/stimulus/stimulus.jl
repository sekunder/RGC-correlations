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

default_CRCNS_dir = ispath("/Users/sekunder/Documents/MATLAB/crcns_ret-1") ? "/Users/sekunder/Documents/MATLAB/crcns_ret-1" : "/data1/homes/abk170/crcns_ret-1"

"""
    framerate(S::AbstractStimulus)

Convenience method returns `1/frametime(S)`. Each concrete subtype of
AbstractStimulus should implement `frametime`.
"""
frame_rate(S::AbstractStimulus) = 1.0 / frametime(S)

function _compute_STRFs(spike_histogram::Matrix{Float64}, stim::AbstractStimulus; kwargs...)
    # TODO implement this. The idea is that frame_vector will return a matrix
    # where each column is one frame of the stimulus, so this function should
    # generically work on any stimulus with frame_vector implemented.
end

################################################################################
#### Include definitions of concrete types
################################################################################
include("CRCNS_Stimulus.jl")

export show, frame_size, frame_time, frame_rate, matrix_form, frame_image,
       CRCNS_Stimulus

end
