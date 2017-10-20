"""
    BernoulliCodeDistribution

Distribution where bits are each independent Bernoulli trials (coin flips)
"""
type BernoulliCodeDistribution <: AbstractBinaryVectorDistribution
    p::Vector{Float64}
    I::Vector{Int}
    metadata::Dict

    function BernoulliCodeDistribution(P::Vector{Float64}, I=collect(1:length(P)); kwargs...)
        new(P, I, Dict(kwargs))
    end
end
(Pr::BernoulliCodeDistribution)(x::Union{Vector{Bool},BitVector}) = length(x) == length(Pr.p) ? prod([x[i] ? Pr.p[i] : (1-Pr.p[i]) for i = 1:length(Pr.p)])

function show(io:IO, P::BernoulliCodeDistribution)
    println(io,"Bernoulli Code Distribution")
    println(io,"N_neurons: $(length(P.p))")
    println(io,"p_i's:     $(P.p)")
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

n_bits(B::BernoulliCodeDistribution) = length(B.p)

function expectation_matrix(B::BernoulliCodeDistribution)
    em = B.p * B.p'
    em[1:(n_bits(B)+1):end] = B.p
    return em
end

entropy(P::BernoulliCodeDistribution) = -sum([(p_i * log(p_i) + (1-p_i) * log(1-p_i)) for p_i in filter(x -> 0 < x < 1, P.p)])

random(B::BernoulliCodeDistribution, n_samples=1) = rand(n_bits(B), n_samples) .< B.p
