
"""
    loglikelihood(X, P)

Computes the log likelihood of `X` given distribution `P`
"""
function loglikelihood(X::Union{Matrix{Bool},BitMatrix}, P::AbstractBinaryVectorDistribution; kwargs...)
    sum_kbn([log(pdf(P,X[:,k])) for k = 1:size(X,2)]) / size(X,2)
end

"""
    loglikelihood(X, Jtilde, grad)

Computes the loglikelihood of the data given the Ising distribution generated by
`Jtilde`. If `grad` is a vector of length > 0, it is modified in-place with the
gradient of the function.

"""
function loglikelihood(X::Union{Matrix{Bool}, BitMatrix}, Jtilde::Vector, grad::Vector=[]; mu_X=(X*X'/size(X,2)), kwargs...)
    (N_neurons, N_samples) = size(X)
    Jtilde = reshape(Jtilde, N_neurons,N_neurons)
    P = IsingDistribution(Jtilde)
    # the way i implemented this before, I just called the pdf function over and
    # over again. So in fact, theoretically the caching I implemented would
    # provide a speed up here. In fact, I should change the _expectation_matrix
    # function so that it invokes get_pdf, so that BCD and ID will cache all the
    # values. Then things should be speedy?
    if length(grad) > 0
        # performing this computation first, if it's necessary, results in all
        # the pdf values being cached.
        # mu_X = X * X' / N_samples
        mu_P = expectation_matrix(P)
        grad[:] = -0.5 * (mu_P - mu_X)[:]
        grad[1:(N_neurons+1):end] = diag(mu_P - mu_X)

    end
    return sum_kbn([log(pdf(P, X[:,k])) for k = 1:N_samples]) / N_samples
end

"""
    MPF_objective(X, J, grad)

Uses function `K` from the MPF paper (or maybe from the MPF sample code for
Matlab). Using this to fit J will typically result in values close to optimal,
but the advantage is that this method does not involve computing the partition
function Z.

"""
function MPF_objective(X::Union{Matrix{Bool}, BitMatrix}, Jtilde::Vector, grad::Vector=[])
    (N_neurons, N_samples) = size(X)
    Jtilde = reshape(Jtilde, N_neurons, N_neurons)
    J = Jtilde - Diagonal(Jtilde)
    theta = diag(Jtilde)
    DeltaX = 2 * X - 1 # this is Δx_l for each codeword x, for each index l
    Kfull = exp((-0.5 * DeltaX .* (J * X) + DeltaX .* theta)/2)
    K = sum_kbn(Kfull[:]) / N_samples
    if length(grad) > 0
        M = zeros(N_neurons, N_neurons)
        M[1:(N_neurons+1):end] = sum(0.5 * Kfull .* DeltaX, 2)
        for p = 1:(N_neurons - 1)
            for q = (p+1):N_neurons
                M[p,q] = M[q,p] = sum([0.5 * -0.5 * (Kfull[p,w] * DeltaX[p,w] * X[q,w] + Kfull[q,w] * DeltaX[q,w] * X[p,w]) for w = 1:N_samples])
            end
        end
        grad[:] = M[:]
    end
    return K
end
