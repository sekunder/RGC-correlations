
"""
    inhomogeneous_poisson_process(lam, dt; kwargs...)

Returns a sequence of spike times produced by the inhomogeneous poisson process
with rate `lambda`. `dt` is 1/sampling rate.

To produce a more reliable spike train, increase `sampling_rate_factor`. If
`sampling_rate_factor = k` then the function first computes `kron(lam, ones(k))`
and uses sampling rate `dt/k` (For instance, with the CRCNS data, `k` between 10
and 20 seems to work).

Will not produce spikes with repeated spike times (this could come up if, e.g.
`lam[t]` is very large)

Options:
 * `sampling_rate_factor::Int = 1` Increases sampling rate

Stopping criteria:
 * `Tmax = length(lam) * dt` Will produces spikes until first spike later than this time
 * `max_real_time` Will stop this process early if it's taking forever (i.e. if `lam` is very low). Default value `0` ignores this option.
 * `max_loops` Same as above, counting number of loops instead of amount of real time. Default `0` to ignore.
 * `max_spikes` Same as above, counting number of spikes produced (e.g. if you want the first-so-many spieks produced by a process.)

These additional kwargs can be passed arrays. The first entry will be modified
in place. This allows additional information to be collected about the process
without requiring the function to return multiple values.
 * `total_time` The total run time of the function
 * `total_loops` The total number of times the process looped
 * `total_spikes` The same as `length(times)` where `times` is the output of this function.

This implementation is based on the first algorithm outlined at
http://freakonometrics.hypotheses.org/724

"""
function inhomogeneous_poisson_process(lam::Vector{Float64}, dt::Float64;
    Tmax::Float64=length(lam) * dt, sampling_rate_factor::Int=1,
    max_real_time::Float64=0.0, max_loops::Int=0, max_spikes::Int=0,
    total_time=[], total_loops=[], total_spikes=[],
    exit_status=Symbol[])
    # based on the first algorithm outlined at http://freakonometrics.hypotheses.org/724
    dt = dt/sampling_rate_factor
    # Lam = cumsum_kbn(kron(lam, ones(sampling_rate_factor)) * dt)
    Lam = cumsum_kbn(kron(lam, fill(dt, sampling_rate_factor)))
    times = Vector{Float64}()
    t_last = 0.0
    s = 0.0; t = 0.0
    real_time = 0.0
    n_loops = 0
    n_spikes = 0
    while t < Tmax
        tic()
        u = 1 - rand()
        s = s - log(u)
        idx = searchsortedlast(Lam, s)
        t = idx * dt
        if Tmax > t > t_last
            # Avoid repeated spikes by comparing to t_last
            push!(times, t)
            t_last = t
            n_spikes += 1
        end
        real_time += toq()
        n_loops += 1
        if idx == length(Lam)
            # if we've hit the end of Lam, then we should quit while we're ahead
            push!(exit_status,:end_of_lambda_reached)
            break
        end
        if 0.0 < max_real_time < real_time
            push!(exit_status,:max_real_time_reached)
            break
        end
        if 0 < max_loops < n_loops
            push!(exit_status,:max_loops_reached)
            break
        end
        if 0 < max_spikes < n_spikes
            push!(exit_status,:max_number_of_spikes)
            break
        end
    end
    if t >= Tmax
        push!(exit_status,:max_spike_time_reached)
    end
    if length(total_time) > 0
        total_time[1] = real_time
    end
    if length(total_loops) > 0
        total_loops[1] = n_loops
    end
    if length(total_spikes) > 0
        total_spikes[1] = length(times)
    end
    return times
end

# MAYBEDO maybe implement a version of this function that returns a proper
# SpikeTrains object. Probably not necessary, but keeping this flag here just in
# case.
