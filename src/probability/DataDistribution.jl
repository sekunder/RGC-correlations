
"""
    DataDistribution

A probability distribution representing the frequency of codewords in a given
sample. The sample must be a `Matrix{Bool}` or a `BitMatrix` with columns
representing codewords. The number of bits must be less than 64

"""
type DataDistribution <: AbstractBinaryVectorDistribution
    P::SparseVector{Float64}
    N::Int
    X::BitMatrix
    I::Vector{Int}
    metadata::Dict

    function DataDistribution(X::Matrix{Bool}, I=collect(1:size(X,1)); kwargs...)
        (N, m) = size(X)
        if N >= 64
            error("Attempting to compute data distribution of data with > 64 bit words")
        end
        Xtilde = sum(X .* [2^i for i = 0:(N - 1)], 1) + 1
        P = sparsevec(Xtilde[:], ones(m), 2^N)
        new(P/m, N, BitMatrix(X), I, Dict(kwargs))
    end
    function DataDistribution(X::BitMatrix, I=collect(1:size(X,1)); kwargs...)
        (N, m) = size(X)
        P = sparsevec([1 + sum(Int, X[:,k].chunks) for k = 1:m], ones(m), 2^N)
        new(P/m, N, X, I, Dict(kwargs))
    end
end
(Pr::DataDistribution)(x::Vector{Bool}) = length(x) == Pr.N ? Pr.P[1 + dot([2^i for i = 0:(Pr.N - 1)], x)] : error("DataDistribution: out of domain error")
(Pr::DataDistribution)(x::BitVector) = length(x) == Pr.N ? Pr.P[1 + sum(Int,x.chunks)] : error("DataDistribution: out of domain error")

function show(io::IO, P::DataDistribution)
    println(io, "Data Distribution")
    println(io, "N_neurons: $(P.N)")
    println(io, "Indices:   $(P.I)")
    println(io, "N_samples: $(size(P.X,2))")
    show_metadata(P)
    # if isempty(P.metadata)
    #     println(io, "No metadata found")
    # else
    #     println(io, "Metadata:")
    #     for (k,v) in P.metadata
    #         println(io, "\t$k : $v")
    #     end
    # end
end

n_bits(DD::DataDistribution) = DD.N

expectation_matrix(DD::DataDistribution) = DD.X * DD.X' / size(DD.X,2)

entropy(DD::DataDistribution) = -sum_kbn([p * log(p) for p in nonzeros(DD.P)])

function random(DD::DataDistribution, n_samples=1)
    cdf = cumsum_kbn(full(DD.P))
    X = falses(n_bits(DD), n_samples)
    r = rand(n_samples)
    for s = 1:n_samples
        k = searchsortedfirst(cdf, r[s])
        X[:,s] = digits(k, 2, n_bits(DD))
    end
    return X
end
