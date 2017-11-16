"""
    BernoulliCodeDistribution

Distribution where bits are each independent Bernoulli trials (coin flips)
"""
type BernoulliCodeDistribution <: AbstractBinaryVectorDistribution
    p::Vector{Float64}
    I::Vector{Int}
    metadata::Dict{Any,Any}
    cache::Dict{Any,Any}

    function BernoulliCodeDistribution(P::Vector{Float64}, I=1:length(P); kwargs...)
        new(P, collect(I), Dict(kwargs), Dict(:pdf=>spzeros(2^length(P))))
    end
end

n_bits(B::BernoulliCodeDistribution) = length(B.p)

function pdf(Pr::BernoulliCodeDistribution, x::Union{Vector{Bool},BitVector})
    if length(x) != n_bits(Pr)
        error("BernoulliCodeDistribution pdf: out of domain error")
    else
        # minor inefficiency: this won't differentiate between a 0.0 value that
        # is stored vs. not stored in the sparse vector. But, this would only
        # come up if one of the p_i = 0.
        # also, might be worth splitting this back into two cases (Vector{Bool}
        # vs. BitVector)
        # idx = 1 + dot([2^i for i = 0:(n_bits(Pr) - 1)], x)
        idx = 1 + _binary_to_int(x)
        if Pr.cache[:pdf][idx] == 0.0
            Pr.cache[:pdf][idx] = prod([x[i] ? Pr.p[i] : (1-Pr.p[i]) for i = 1:length(Pr.p)])
        end
        return Pr.cache[:pdf][idx]
        # if Pr.cache[:pdf][idx] > 0.0
        #     return Pr.cache[:pdf][idx]
        # else
        #     Pr.cache[:pdf][idx] = prod([x[i] ? Pr.p[i] : (1-Pr.p[i]) for i = 1:length(Pr.p)])
        #     return Pr.cache[:pdf][idx]
        # end
    end
    # length(x) == length(Pr.p) ? prod([x[i] ? Pr.p[i] : (1-Pr.p[i]) for i = 1:length(Pr.p)]) : error("BernoulliCodeDistribution: out of domain error")
end



function show(io::IO, P::BernoulliCodeDistribution)
    println(io,"Bernoulli Code Distribution")
    println(io,"N_neurons: $(n_bits(P))")
    println(io,"p_i's:     $(P.p)")
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

function expectation_matrix(B::BernoulliCodeDistribution)
    em = B.p * B.p'
    em[1:(n_bits(B)+1):end] = B.p
    return em
end

# I'm going to be lazy here and allow this value to be computed each time. The
# parentheses around the second assignment expression are just for human
# readability.
entropy(P::BernoulliCodeDistribution) = (P.metadata[:entropy] = -sum([(p_i * log(p_i) + (1-p_i) * log(1-p_i)) for p_i in filter(x -> 0 < x < 1, P.p)]))
entropy2(P::BernoulliCodeDistribution) = (P.metadata[:entropy] = -sum([(p_i * log2(p_i) + (1-p_i) * log2(1-p_i)) for p_i in filter(x -> 0 < x < 1, P.p)]))

random(B::BernoulliCodeDistribution, n_samples=1) = rand(n_bits(B), n_samples) .< B.p
