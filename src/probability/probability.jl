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
#### Miscellaneous functions
################################################################################

################################################################################
#### Include definitions of concrete distributions
################################################################################

include("DataDistribution.jl")

export show, n_bits, random, expectation_matrix, entropy,
       DataDistribution

end
