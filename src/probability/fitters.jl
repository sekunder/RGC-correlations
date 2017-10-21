
"""
    first_order_model(X, I; kwargs...)

Returns a `BernoulliCodeDistribution` with `p_i` equal to the expected value of
bit i.

"""
function first_order_model(X::Union{Matrix{Bool},BitMatrix}, I=1:size(X,1); kwargs...)
    mu = sum(X[I,:],2) / length(I)
    return BernoulliCodeDistribution(mu, I; kwargs...)
end

"""
    second_order_model(X, I; kwargs...)

Returns an `IsingDistribution` which is fit to the pairwise correlations in `X`.
"""
function second_order_model(X::Union{Matrix{Bool},BitMatrix}, I=1:size(X,1); kwargs...)
    # TODO fit the Ising model, using either loglikelihood or MPF
end

"""
    data_model(X, I; kwargs...)

Returns a `DataDistribution` based on the specified data and indices.
"""
function data_model(X::Union{Matrix{Bool},BitMatrix}, I=1:size(X,1); kwargs...)
    # TODO this should be pretty straightforward
end
