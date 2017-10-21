
using NLopt

"""
    first_order_model(X, I; kwargs...)

Returns a `BernoulliCodeDistribution` with `p_i` equal to the expected value of
bit i.

"""
function first_order_model(X::Union{Matrix{Bool},BitMatrix}, I=1:size(X,1); kwargs...)
    mu = sum(X[I,:],2) / length(I)
    return BernoulliCodeDistribution(mu, I; autocmment="first_order_model", kwargs...)
end

"""
    second_order_model(X, I; kwargs...)

Returns an `IsingDistribution` which is fit to the pairwise correlations in `X`.
"""
function second_order_model(X::Union{Matrix{Bool},BitMatrix}, I=1:size(X,1); kwargs...)
    Xselected = X[I,:]
    N_neurons,N_samples = size(Xselected)

    Jseed = rand(N_neurons,N_neurons); Jseed = (Jseed + Jseed') / (2 * N_neurons)

    # let's try only writing one method, since the only difference is the
    # function and max/min.
    opt_Ising = Opt(:LD_LBFGS, N_neurons^2)
    if N_neurons <= ISING_METHOD_THRESHOLD
        #TODO current implementation will, say, compute mu_X repeatedly. Change loglikelihood to accept mu as an optional parameter.
        F_X(J, g) = loglikelihood(Xselected, J, g)
        max_objective!(opt_Ising, F_X)
        fun = "loglikelihood"
    else
        F_X(J, g) = MPF_objective(Xselected, J, g)
        min_objective!(opt_Ising, F_X)
        fun = "MPF"
    end
    ftol_rel!(opt_Ising, 1e-10)
    ftol_abs!(opt_Ising, 1e-10)
    xtol_rel!(opt_Ising, 1e-10)
    xtol_abs!(opt_Ising, 1e-10)
    vector_storage!(opt_Ising, N_neurons)
    (optVal, J_opt, optReturn) = optimize(opt_Ising, Jseed[:])
    J_opt = reshape(J_opt, N_neurons, N_neurons)
    return IsingDistribution(J_opt, I; autocomment="second_order_model (using $fun)", opt_val=optVal, opt_ret=optReturn)
end

"""
    data_model(X, I; kwargs...)

Returns a `DataDistribution` based on the specified data and indices.
"""
function data_model(X::Union{Matrix{Bool},BitMatrix}, I=1:size(X,1); kwargs...)
    return DataDistribution(X, I; autocomment="data_model", kwargs...)
end
