
using NLopt
include("../util/optimizers.jl")

"""
    first_order_model(X, I; kwargs...)

Returns a `BernoulliCodeDistribution` with `p_i` equal to the expected value of
bit i.

"""
function first_order_model(X::Union{Matrix{Bool},BitMatrix}, I=1:size(X,1); kwargs...)
    mu = sum(X[I,:],2) / size(X,2)
    return BernoulliCodeDistribution(mu[:], I; autocomment="first_order_model", kwargs...)
end

"""
    second_order_model(X, I; kwargs...)

Returns an `IsingDistribution` which is fit to the pairwise correlations in `X`.
"""
function second_order_model(X::Union{Matrix{Bool},BitMatrix}, I=1:size(X,1); verbose=false, kwargs...)
    dkwargs = Dict(kwargs)
    if get(dkwargs, :algorithm, :naive) == :naive
        delete!(dkwargs, :algorithm)
        return _Naive_second_order_model(X, I; verbose=verbose, dkwargs...)
    else
        return _NLopt_second_order_model(X, I; verbose=verbose, dkwargs...)
    end
end
function _Naive_second_order_model(X::Union{Matrix{Bool},BitMatrix}, I=1:size(X,1); verbose=false, kwargs...)
    dkwargs = Dict(kwargs)
    println("X is a $(typeof(X)), size(X) = $(size(X))") #TODO debugging
    Xselected = X[I,:]
    N_neurons,N_samples = size(Xselected)
    println("Xselected is a $(typeof(Xselected)), size(Xselected) = $(size(Xselected))") #TODO debugging

    Jseed = rand(N_neurons,N_neurons); Jseed = (Jseed + Jseed') / (2 * N_neurons)
    mu = Xselected * Xselected' / N_samples
    println("mu is a $(typeof(mu)), size(mu) = $(size(mu))") #TODO debugging
    L_X(J,g) = loglikelihood(Xselected, J, g; mu_X = mu) #TODO debugging
    fun = "loglikelihood"
    objective = :max
    if N_neurons > ISING_METHOD_THRESHOLD || pop!(dkwargs, :force_MPF, false)
        K_X(J,g) = MPF_objective(Xselected, J, g) #TODO debugging
        fun = "MPF"
        objective = :min
    end

    if verbose
        println("second_order_model[gradient_optimizer]: setting $objective objective $fun")
    end
    verbosity = verbose + pop!(dkwargs, :more_verbose, false)
    # (F_opt, J_opt, stop) = gradient_optimizer(F_X, Jseed[:]; objective=objective, verbose=verbosity, dkwargs...)
    # final_val = F_X(J_opt, [])  #TODO debugging
    if fun == "loglikelihood"
        (F_opt, J_opt, stop) = gradient_optimizer(L_X, Jseed[:]; objective=objective, verbose=verbosity, dkwargs...)
    else
        (F_opt, J_opt, stop) = gradient_optimizer(K_X, Jseed[:]; objective=objective, verbose=verbosity, dkwargs...)
    end
    J_opt = reshape(J_opt, N_neurons, N_neurons)
    return IsingDistribution(J_opt, I; autocomment="second_order_model[gradient_optimizer|$fun]", opt_val=F_opt, opt_ret=stop, dkwargs...)
    # return (F_opt, J_opt, stop, Jseed, mu, (fun == "loglikelihood" ? L_X : K_X))
end
function _NLopt_second_order_model(X::Union{Matrix{Bool},BitMatrix}, I=1:size(X,1); verbose=false, kwargs...)
    dkwargs = Dict(kwargs)
    Xselected = X[I,:]
    N_neurons,N_samples = size(Xselected)

    Jseed = rand(N_neurons,N_neurons); Jseed = (Jseed + Jseed') / (2 * N_neurons)

    # for some reaosn I'm getting "F_X not defined". So, I'm going to try
    # definind F_X as loglikelihood no matter what, then override F_X with MPF
    # if necessary. I'll set a flag for which function I used, and then set
    # max_objective or min_objective as necessary.
    # Update: turns out maybe Julia was setting F_X no matter whether it
    # actually entered the branch? Or something.
    mu = Xselected * Xselected' / N_samples
    L_X(J,g) = loglikelihood(Xselected, J, g; mu_X=mu)
    fun = "loglikelihood"
    if N_neurons > ISING_METHOD_THRESHOLD || pop!(dkwargs, :force_MPF, false)
        K_X(J,g) = MPF_objective(Xselected, J, g)
        fun = "MPF"
    end

    # let's try only writing one method, since the only difference is the
    # function and max/min.
    alg = pop!(dkwargs,:algorithm,:LD_LBFGS)
    opt_Ising = Opt(alg, N_neurons^2)

    if fun == "loglikelihood"
        if verbose
            println("second_order_model[NLopt/$alg]: setting max objective function $fun")
        end
        max_objective!(opt_Ising, L_X)
    else
        if verbose
            println("second_order_model[NLopt/$alg]: setting min objective function $fun")
        end
        min_objective!(opt_Ising, K_X)
    end
    if haskey(dkwargs, :ftol_rel)
        ftol_rel!(opt_Ising, pop!(dkwargs,:ftol_rel))
    end
    if haskey(dkwargs, :ftol_abs)
        ftol_abs!(opt_Ising, pop!(dkwargs,:ftol_abs))
    end
    if haskey(dkwargs, :xtol_rel)
        xtol_rel!(opt_Ising, pop!(dkwargs,:xtol_rel))
    end
    if haskey(dkwargs, :xtol_abs)
        xtol_abs!(opt_Ising, pop!(dkwargs,:xtol_abs))
    end
    if haskey(dkwargs, :vector_storage)
        vector_storage!(opt_Ising, pop!(dkwargs,:vector_storage))
    end
    if haskey(dkwargs, :maxeval)
        maxeval!(opt_Ising, pop!(dkwargs,:maxeval))
    end
    if verbose
        println("second_order_model[NLopt/$alg]: running optimization")
        println("\talgorithm: $(algorithm(opt_Ising))")
        println("\tftol (rel/abs): $(ftol_rel(opt_Ising)) / $(ftol_abs(opt_Ising))")
        println("\txtol (rel/|abs|): $(xtol_rel(opt_Ising)) / $(norm(xtol_abs(opt_Ising)))")
        println("\tvect. stor.: $(vector_storage(opt_Ising))")
        println("\tmax evals: $(maxeval(opt_Ising))")
    end
    (optVal, J_opt, optReturn) = optimize(opt_Ising, Jseed[:])
    # final_val = F_X(J_opt, [])
    J_opt = reshape(J_opt, N_neurons, N_neurons)
    return IsingDistribution(J_opt, I; autocomment="second_order_model[NLopt/$alg|$fun]", opt_val=optVal, opt_ret=optReturn, dkwargs...)
end

"""
    data_model(X, I; kwargs...)

Returns a `DataDistribution` based on the specified data and indices.
"""
function data_model(X::Union{Matrix{Bool},BitMatrix}, I=1:size(X,1); kwargs...)
    return DataDistribution(X, I; autocomment="data_model", kwargs...)
end
