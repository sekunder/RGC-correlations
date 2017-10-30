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
    ISING_METHOD_THRESHOLD

If N_neurons is greater than this number, switch from exact sampling to Gibbs
sampling and use MPF instead of LogLikelihood to for the parameters.
"""
const ISING_METHOD_THRESHOLD = 20
const METH_THRESH_WARN_KEY = 1

function show_metadata(io::IO, P::AbstractBinaryVectorDistribution)
    if isempty(P.metadata)
        println(io, "No metadata found")
    else
        println(io, "Metadata:")
        for (k,v) in P.metadata
            println(io, "\t$k : $v")
        end
    end
end

"""
    binary_to_int(x::Union{BitVector, Vector{Bool}})

Returns the integer that `x` represents. Will use `x.chunks` for `BitVector`s
and the "dot product method" for `Vector{Bool}`s. Note that this will really
mess you up if `length(x) > 63`

"""
_binary_to_int(x::BitVector) = Int(x.chunks[1])
_binary_to_int(x::Vector{Bool}) = dot([2^i for i = 0:(length(x) - 1)], x)

"""
    get_pdf(P), get_cdf(P)

Internally-used function to return a vector of the PDF of the distribution.
Naive implementation just calls pdf(P,x) over all possible x. Note that for this
internal method to work correctly, pdf(P,x) should cache the value in
P.cache[:pdf] as a side effect.
"""
function _get_pdf(P::AbstractBinaryVectorDistribution)
    if get(P.metadata, :pdf_computed, false)
        return full(P.cache[:pdf])
    else
        P.metadata[:pdf_computed] = true
        return [pdf(P, digits(Bool,x,2,n_bits(P))) for x in 0:(2^n_bits(P) - 1)]
    end
end

function _get_cdf(P::AbstractBinaryVectorDistribution)
    if haskey(P.cache, :cdf)
        cdf = P.cache[:cdf]
    else
        cdf = cumsum_kbn(get_pdf(P))
        P.cache[:cdf] = cdf
    end
    return cdf
end

"""
    random_exact(P, n_sample)

Internally-used function for using a naive sampling method to draw random
vectors from a distribution. Relies on an implementation of `get_cdf` for the
distribution (so, don't use with distributions on more than ~20 neurons)

"""
function _random_exact(P::AbstractBinaryVectorDistribution, n_samples::Int=1)
    X = falses(n_bits(P), n_samples)
    r = rand(n_samples)
    cdf = get_cdf(P)
    for s = 1:n_samples
        k = searchsortedfirst(cdf, r[s])
        X[:,s] = digits(Bool, k - 1, 2, n_bits(P))
    end
    return X
end

"""
    expectation_matrix(P)

Returns a matrix with the expected values of bits and pairs of bits. Default
behavior is to call get_pdf(P) and then loop through those values.

"""
function _expectation_matrix(P::AbstractBinaryVectorDistribution)
    em = zeros(n_bits(P), n_bits(P))
    p = get_pdf(P)
    for k = 0:(2^n_bits(P) - 1)
        x = digits(Bool, k, 2, n_bits(P))
        em += p[k+1] * x * x'
    end
    return em
end

_entropy(P::AbstractBinaryVectorDistribution) = -sum_kbn([p * log(p) for p in filter(x -> 0 < x < 1, get_pdf(P))])

################################################################################
#### Include definitions of concrete distributions, other functions
################################################################################
include("DataDistribution.jl")
include("BernoulliCodeDistribution.jl")
include("IsingDistribution.jl")

include("fitters.jl")
include("loglikelihood.jl")

export AbstractBinaryVectorDistribution,
       show, n_bits, random, expectation_matrix, entropy, pdf, get_pdf, get_cdf,
       DataDistribution, BernoulliCodeDistribution, IsingDistribution,
       first_order_model, second_order_model, data_model,
       loglikelihood, MPF_objective

end
