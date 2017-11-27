include("../util/optimizers.jl")
"""
    STRF_response(STRF, stimulus; kwargs...)

Computes the "raw" response rate of a neuron with spatio-temporal field given by
`STRF` in response to `stimulus`. Returns a pair `(r,tau)` where `r` is a vector
of rate values and `tau` is the temporal resolution (i.e. the width in seconds
of each time bin represented by `r`). For now, I assume the STRF and the
stimulus use the same frame rate, meaning `tau = frame_time(STRF)`. I'm leaving
the door open to needing to fiddle with this in the future, though (i.e. if
we're comparing an STRF computed at 30Hz to a stimulus recorded at 60Hz, that
sort of thing.)

The "raw" response rate is given by convolving `STRF` and `stimulus` as
functions of time.

"""
function STRF_response(STRF::AbstractStimulus, stimulus::AbstractStimulus;
    flip_STRF_time=false, flip_stimulus_time=false,
    kwargs...)

    M_strf = matrix_form(STRF); M_stim = matrix_form(stimulus)
    if flip_STRF_time
        M_strf = flipdim(M_strf, 2)
    end
    if flip_stimulus_time
        M_stim = flipdim(M_stim, 2)
    end
    dx = prod(STRF.d .* STRF.mm_per_px)
    N_bits_per_frame, N_frames = size(M_stim)
    r = zeros(N_frames)
    for bit = 1:N_bits_per_frame
        r[:] += conv(M_strf[bit,:], M_stim[bit,:])[1:length(r)]
    end
    return r * dx * frame_time(STRF), frame_time(STRF)
end

"""
    scale_response(r, n, phi, Q; objective=:min, kwargs...)

Finds parameters for the function `phi` to minimize penalty function `Q`.
Returns the scaled output, the optimal parameters, and the value of `Q`.

Mathematically, I'm thinking of `ϕ(x; θ)` as the transfer function which I'm
applying to `r`, and `Q` as the penalty function for how close `ϕ(r(t); θ)` is
to `n(t)`.

Programmatically, this means:
 * `r::Vector{Float64}` is probably a rate computed from `STRF_response`
 * `n::Vector{Float64}` is probably a single neuron's spike histogram
 * `phi(x, theta)` should accept vectors `x` and `theta` and return a vector the same length as `x`
 * `Q(L,n)` should accept two vectors and return a single number

Some keyword arguments:
 * `objective=:min` whether to maximizer or minimize `Q`
 * `d` is the dimensionality of the search space (if using brute force search)
 * `theta0` is an initial guess for the value of the parameters (if using some other algorithm)
 * `algorithm=:brute_force_search` how to search for the optimum

Currently, only supports brute force search, but I'm leaving the door open to
expand this if it seems necessary. (If we're being honest this function is
mostly unnecessary anyway since all it's doing is acting as a wrapper to...
another function. But it'll make my scripts readable so whatever.)

See the `brute_force_optimizer` documentation for some necessary kwargs.

"""
function scale_response(r::Vector{Float64}, n::Vector{Float64}, phi, Q;
    algorithm=:brute_force_search, kwargs...)

    # dkwargs = Dict(kwargs)
    if algorithm != :brute_force_search
        #MAYBEDO expand to other options
        error("scale_response: No, seriously, the only algorithm option is :brute_force_search (got $algorithm)")
    end
    # d = pop!(dkwargs, :d, 3)
    obj_fun(x::Vector) = Q(phi(r, x), n)
    Q_opt, theta_opt, criteria = brute_force_optimizer3(obj_fun, 3; kwargs...)
    return phi(r, theta_opt), theta_opt, Q_opt
end
