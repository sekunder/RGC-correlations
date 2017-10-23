"""
Spikes

Includes types and methods for handling spike trains.
"""
module Spikes

import Base.show

################################################################################
#### Include appropriate files
################################################################################
include("SpikeTrains.jl")
include("poissonproc.jl")

################################################################################
#### Miscellaneous functions and constants
################################################################################

function show_metadata(io::IO, ST::SpikeTrains)
    if isempty(ST.metadata)
        println(io, "No metadata found")
    else
        println(io, "Metadata:")
        for (k,v) in ST.metadata
            println(io, "\t$k : $v")
        end
    end
end

export show,
       SpikeTrains, histogram, raster,
       inhomogeneous_poisson_process

end
