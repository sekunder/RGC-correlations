
"""
    DataDistribution(X, I; kwargs...)

A probability distribution representing the frequency of codewords in a given
sample. The sample must be a `Matrix{Bool}` or a `BitMatrix` with columns
representing codewords. The number of bits = `length(I)` must be less than 64.

First sets `X = X[I,:]` then proceeds

"""
type DataDistribution <: AbstractBinaryVectorDistribution
    P::SparseVector{Float64}
    N::Int
    X::BitMatrix
    I::Vector{Int}
    metadata::Dict{Any,Any}
    cache::Dict{Any,Any}

    function DataDistribution(X::Union{Matrix{Bool},BitMatrix}, I=1:size(X,1); kwargs...)
        X_subset = X[I,:]
        (N,m) = size(X_subset)
        if N >= 64
            error("Attempting to compute data distribution of data with >= 64 bit words")
        end
        P = sparsevec(_compute_Xtilde(X_subset)[:], ones(m), 2^N)
        new(P/m, N, BitMatrix(X_subset), collect(I), Dict(kwargs), Dict())
    end
end
_compute_Xtilde(X::Matrix{Bool}) = 1 .+ sum(X .* [2^i for i = 0:(size(X,1)-1)], 1)
_compute_Xtilde(X::BitMatrix) = [1 + Int(X[:,k].chunks[1]) for k = 1:size(X,2)]


n_bits(DD::DataDistribution) = DD.N

# pdf(Pr::DataDistribution, x::Vector{Bool}) = length(x) == n_bits(Pr) ? Pr.P[1 + dot([2^i for i = 0:(Pr.N - 1)], x)] : error("DataDistribution pdf: out of domain error")
# pdf(Pr::DataDistribution, x::BitVector) = length(x) == n_bits(Pr) ? Pr.P[1 + Int(x.chunks[1])] : error("DataDistribution pdf: out of domain error")
pdf(Pr::DataDistribution, x::Union{Vector{Bool},BitVector}) = length(x) == n_bits(Pr) ? Pr.P[1 + _binary_to_int(x)] : error("DataDistribution pdf: out of domain error")

get_pdf(DD::DataDistribution) = full(DD.P)


function show(io::IO, P::DataDistribution)
    println(io, "Data Distribution")
    println(io, "N_neurons: $(n_bits(P))")
    println(io, "Indices:   $(P.I)")
    println(io, "N_samples: $(size(P.X,2))")
    show_metadata(io, P)
    # if isempty(P.metadata)
    #     println(io, "No metadata found")
    # else
    #     println(io, "Metadata:")
    #     for (k,v) in P.metadata
    #         println(io, "\t$k : $v")
    #     end
    # end
end

expectation_matrix(DD::DataDistribution) = DD.X * DD.X' / size(DD.X,2)

entropy(DD::DataDistribution) = -sum_kbn([p * log(p) for p in nonzeros(DD.P)])
entropy2(DD::DataDistribution) = -sum_kbn([p * log2(p) for p in nonzeros(DD.P)])

# function random(DD::DataDistribution, n_samples=1)
#     X = falses(n_bits(DD), n_samples)
#     r = rand(n_samples)
#     for s = 1:n_samples
#         k = searchsortedfirst(get_cdf(DD), r[s])
#         X[:,s] = digits(k, 2, n_bits(DD))
#     end
#     return X
# end
random(DD::DataDistribution, n_samples::Int=1) = _random_exact(DD, n_samples)
