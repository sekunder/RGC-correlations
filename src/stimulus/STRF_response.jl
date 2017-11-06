
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
    dx = prod(STRF.d) * STRF.mm_per_px
    N_bits_per_frame, N_frames = size(M_stim)
    r = zeros(N_frames)
    for bit = 1:N_bits_per_frame
        ans[:] = conv(M_strf[bit,:], M_stim[bit,:])[1:length(r)]
    end
    return r * dx * frame_time(STRF), frame_time(STRF)
end

# TODO why don't I just write the generic brute force search optimizer
# TODO then I can define some penalty functions here
# """
#     scale_response(r, n, phi, Q, theta0; objective=:min, lb=0.0, ub=1.0, kwargs...)
#
# Finds parameters for the function `phi` to minimize penalty function `Q`,
# """
# function scale_response(r::Vector{Float64}, n::Vector{Float64}, phi, Q, theta0::Vector{Float64};
#     objective=:max
#     lb::Union{Float64,Vector{Float64}}=0.0, ub::Union{Float64,Vector{Float64}}=1.0,
#     kwargs...)
# end
