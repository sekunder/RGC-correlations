"""
Probability

Includes types and methods for representing probability distributions on the set of binary vectors.
"""
module Probability

import Base.show

################################################################################
#### Abstract type to use as placeholder/to get all docstrings in one place
################################################################################
"""
    AbstractBinaryVectorDistribution

Abstract type used as placeholder.
"""
abstract AbstractBinaryVectorDistribution

################################################################################
#### Miscellaneous functions and constants
################################################################################

"""
    ISING_SAMPLE_METHOD_THRESHOLD

If N_neurons is greater than this number, switch from exact sampling to Gibbs
sampling and use MPF instead of LogLikelihood to for the parameters.
"""
const ISING_SAMPLE_METHOD_THRESHOLD = 20

function show_metadata(io, P<:AbstractBinaryVectorDistribution)
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
#### Include definitions of concrete distributions
################################################################################

include("DataDistribution.jl")
include("BernoulliCodeDistribution.jl")
include("IsingDistribution.jl")

export show, n_bits, random, expectation_matrix, entropy,
       DataDistribution, BernoulliCodeDistribution, IsingDistribution

end
