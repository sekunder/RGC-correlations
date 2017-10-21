
"""
    inhomogeneous_poisson_process(lam, dt, tmax = length(lam) * dt)

    Returns a sequence of spike times produced by the inhomogeneous poisson
    process with rate lambda. dt is 1/sampling rate. To produce a more reliable
    spike train, increase the sampling rate. This can be achieved by
    `NaiveInhomogPoisson(kron(lam, ones(k)), dt/k)` (For instance, with the
    CRCNS data, `k` between 10 and 20 seems to work)

    Will not produce spikes with repeated spike times (this could come up if,
    e.g. `lam[t]` is close to `0`)

"""
function inhomogeneous_poisson_process(lam::Vector{Float64}, dt::Float64, Tmax::Float64=length(lam) * dt)
  # based on the first algorithm outlined at http://freakonometrics.hypotheses.org/724
  Lam = cumsum_kbn(lam * dt)
  times = Vector{Float64}()
  t_last = 0.0
  s = 0.0; t = 0.0
  while t < Tmax
    u = (1 - rand(1))[1]
    s = s - log(u)
    idx = searchsortedlast(Lam, s)
    # t = (idx + (s - Lam[idx])/lam[idx]) * dt
    t = idx * dt
    # t = lam[idx] == 0 ? idx * dt : (idx + (s - Lam[idx])/lam[idx]) * dt
    # if idx > 1 && lam[idx] != 0
    #   t += (s - Lam[idx])/lam[idx]
    # end
    if t > t_last
      push!(times, t)
      t_last = t
    end
  end
  return times
end

# TODO maybe implement a version of this function that returns a proper
# SpikeTrains object. Probably not necessary, but keeping this flag here just in
# case.
