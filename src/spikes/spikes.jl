"""
Spikes

Includes types and methods for handling spike trains.
"""
module Spikes

import Base.show

################################################################################
#### Miscellaneous functions and constants
################################################################################

function show_metadata(io::IO, P::AbstractBinaryVectorDistribution)
    if isempty(P.metadata)
        println(io, "No metadata found")
    else
        println(io, "Metadata:")
        for (k,v) in P.metadata
            # TODO It might be useful to store things like the pdf/cdf in metadata
            # so, think about how to handle that
            println(io, "\t$k : $v")
        end
    end
end

################################################################################
#### Include appropriate files
################################################################################
include("SpikeTrains.jl")
include("poissonproc.jl")

export show,
       SpikeTrains, histogram, raster,
       inhomogeneous_poisson_process

end
