"""
Spikes

Includes types and methods for handling spike trains.
"""
module Spikes

import Base.show
include("../util/metadata.jl")

################################################################################
#### Include appropriate files
################################################################################
include("SpikeTrains.jl")
include("poissonproc.jl")
include("CRCNS/CRCNS.jl")

################################################################################
#### Miscellaneous functions and constants
################################################################################

export show,
       SpikeTrains, histogram, raster, n_cells,
       inhomogeneous_poisson_process

end
