"""
    IsingDistribution

Represents the Ising distribution `P(x) = 1/Z exp(-1/2 x' J x + x' th)`
"""
type IsingDistribution <: AbstractBinaryVectorDistribution
    J::Matrix{Float64}
    theta::Vector{Float64}
    I::Vector{Int}
    metadata::Dict

    function IsingDistribution(J::Matrix{Float64}, theta::Vector{Float64}, I=collect(1:length(theta)); kwargs...)
        if size(J,1) != size(J,2)
            error("Ising Distribution: J must be square")
        end
        if size(J,1) != length(theta)
            error("Ising Distribution: J and theta must have compatible size")
        end
        Jsym = (J + J') / 2
        Jsym = Jsym - Diagonal(Jsym)
        new(Jsym, theta, I, Dict(kwargs))
    end
    function IsingDistribution(Jtilde::Matrix{Float64}, I=collect(1:size(Jtilde,1)); kwargs...)
        if size(J,1) != size(J,2)
            error("Ising Distribution: J must be square")
        end
        theta = diag(Jtilde)
        Jnodiag = Jtilde - Diagonal(Jtilde)
        Jsym = (Jnodiag + Jnodiag') / 2
        new(Jsym, theta, I, metadata)
    end
end
# TODO implement probability calculation

n_bits(P::IsingDistribution) = length(P.theta)

function show(io::IO, P::IsingDistribution)
    println(io, "Ising Distribution")
    println(io, "N_neurons: $(n_bits(P))")
    println(io, "Indices:   $(P.I)")
    show_metadata(P)
end

function expectation_matrix(ID::IsingDistribution)
    # TODO if N_neurons is small enough, we can do this computation
end

function entropy(ID::IsingDistribution)
    # TODO if N_neurons is small enough, we can do this computation
end

function random(ID::IsingDistribution, n_samples=1)
    # TODO if N_neurons is small enough, use exact sampling. otherwise, gibbs.
end
function _random_exact(ID::IsingDistribution, n_samples=1)
    # TODO compute CDF, use searchsortedfirst
end
function _random_gibbs(ID::IsingDistribution, n_samples=1)
    # TODO implement gibbs sampling. Use ID.metadata to decide burnin, independent steps, etc.
end
