"""
    IsingDistribution

Represents the Ising distribution `P(x) = 1/Z exp(1/2 x' J x - x' th)`
"""
type IsingDistribution <: AbstractBinaryVectorDistribution
    J::Matrix{Float64}
    theta::Vector{Float64}
    I::Vector{Int}
    metadata::Dict{Any,Any}
    cache::Dict{Any,Any}

    function IsingDistribution(J::Matrix{Float64}, theta::Vector{Float64}, I=1:length(theta); kwargs...)
        if size(J,1) != size(J,2)
            error("Ising Distribution: J must be square")
        end
        if size(J,1) != length(theta)
            error("Ising Distribution: J and theta must have compatible size")
        end
        Jsym = (J + J') / 2
        Jsym = Jsym - Diagonal(Jsym)
        new(Jsym, theta, collect(I), Dict(kwargs), Dict(:pdf=>spzeros(2^length(theta))))
    end
    function IsingDistribution(Jtilde::Matrix{Float64}, I=collect(1:size(Jtilde,1)); kwargs...)
        theta = diag(Jtilde)
        Jnodiag = Jtilde - Diagonal(Jtilde)
        IsingDistribution(Jnodiag, theta, I; kwargs...)
        # if size(Jtilde,1) != size(Jtilde,2)
        #     error("Ising Distribution: J must be square")
        # end
        # theta = diag(Jtilde)
        # Jnodiag = Jtilde - Diagonal(Jtilde)
        # Jsym = (Jnodiag + Jnodiag') / 2
        # new(Jsym, theta, collect(I), Dict(kwargs), Dict(:pdf=>spzeros(2^length(theta))))
    end
end

################################################################################
#### Miscellaneous computations/internal functions
################################################################################
_E_Ising(ID::IsingDistribution, x::Union{Vector{Bool}, BitVector}) = dot(x, -0.5*ID.J*x + ID.theta)
_get_energies(ID::IsingDistribution) = [_E_Ising(ID, digits(Bool,x,2,n_bits(ID))) for x in 0:(2^n_bits(ID) - 1)]
function _get_Z(ID::IsingDistribution)
    if !haskey(ID.cache, :Z)
        if n_bits(ID) > ISING_METHOD_THRESHOLD
            warn("Computing Z for Ising PDF with $(n_bits(ID)) neurons. This might take a while. (Warning only displays once)", once=true, key=METH_THRESH_WARN_KEY)
        end
        # ID.cache[:Z] = sum_kbn([exp(-_E_Ising(ID, digits(Bool, k, 2, n_bits(ID)))) for k = 0:(2^n_bits(ID) - 1)])
        ID.cache[:Z] = sum_kbn(exp.(-_get_energies(ID)))
    end
    return ID.cache[:Z]
end

################################################################################
#### Implementation of methods for distributions
################################################################################
function entropy(P::IsingDistribution)
    if !haskey(P.metadata, :entropy)
        energies = _get_energies(P)
        # P.metadata[:entropy] = log(_get_Z(P)) + sum_kbn([exp(-_E_Ising(P,digits(Bool,2,x,n_bits(P)))) for x in 0:(2^n_bis(P) - 1)])
        P.metadata[:entropy] = log(_get_Z(P)) + sum_kbn(exp.(-energies) .* energies) / _get_Z(P)
    end
    return P.metadata[:entropy]
end
function entropy2(P::IsingDistribution)
    if !haskey(P.metadata, :entropy2)
        energies = _get_energies(P)
        # P.metadata[:entropy] = log(_get_Z(P)) + sum_kbn([exp(-_E_Ising(P,digits(Bool,2,x,n_bits(P)))) for x in 0:(2^n_bis(P) - 1)])
        P.metadata[:entropy2] = log2(_get_Z(P)) + log2(e) * sum_kbn(exp.(-energies) .* energies) / _get_Z(P)
    end
    return P.metadata[:entropy2]
end

n_bits(P::IsingDistribution) = length(P.theta)

function pdf(ID::IsingDistribution, x::Union{Vector{Bool},BitVector})
    if length(x) != n_bits(ID)
        error("IsingDistribution pdf: out of domain error")
    else
        # if n_bits(ID) > ISING_METHOD_THRESHOLD
        #     warn("Calling Ising PDF with $(n_bits(ID)) neurons. This might take a while. (Warning only displays once)", once=true, key=METH_THRESH_WARN_KEY)
        # end
        # idx = 1 + dot([2^i for i = 0:(n_bits(ID) - 1)], x)
        idx = 1 + _binary_to_int(x)
        if ID.cache[:pdf][idx] == 0.0
            ID.cache[:pdf][idx] = exp(-_E_Ising(ID,x)) / _get_Z(ID)
        end
        return ID.cache[:pdf][idx]
        # if ID.cache[:pdf][idx] > 0.0
        #     p = ID.cache[:pdf][idx]
        # else
        #     p = exp(-_E_Ising(ID, x)) / _get_Z(ID)
        #     ID.cache[:pdf][idx] = p
        # end
        # return p
    end
end


function show(io::IO, P::IsingDistribution)
    println(io, "Ising Distribution")
    println(io, "N_neurons: $(n_bits(P))")
    println(io, "Indices:   $(P.I)")
    show_metadata(io, P)
end


function random(ID::IsingDistribution, n_samples=1, force_gibbs=false)
    if n_bits(ID) <= ISING_METHOD_THRESHOLD && !force_gibbs
        return _random_exact(ID, n_samples)
    else
        return _random_gibbs(ID, n_samples)
    end
end
function _random_gibbs(ID::IsingDistribution, n_samples=1)
    burnin = _get_burnin(ID)
    ind_steps = _get_ind_steps(ID)
    N_steps = burnin + (n_samples - 1) * ind_steps

    d_i = rand(1:n_bits(ID), N_steps) # sequence of bits to update
    u_r = rand(N_steps) # sequence of uniform random samples from [0,1]
    X_out = falses(n_bits(ID), n_samples)
    x = rand(Bool, n_bits(ID)) # the current state

    sample = 1
    next_sample_step = burnin

    for step = 1:N_steps
        # compute the conditional probability of bit d_i[step] given the other bits
        E_act = dot(x, -0.5*ID.J[d_i[step],:]) + ID.theta[d_i[step]]
        P_act = 1 / (1 + exp(-E_act))
        x[d_i[step]] = P_act < u_r[step]
        if step == next_sample_step
            X_out[:,sample] = x
            sample += 1
            next_sample_step += ind_steps
        end
    end
    return X_out
end
_get_burnin(ID::IsingDistribution) = get!(ID.metadata, :burnin, 100*n_bits(ID))
_get_ind_steps(ID::IsingDistribution) = get!(ID.metadata, :ind_steps, 10*n_bits(ID))
