"""
    SpikeTrains(TT, I;
        rec_start, rec_finish, stim_start, stim_finish,
        kwargs...)

A wrapper object for spike trains. `TT` are the spike times, `I` are the
indices/identities of the neurons.

"""
type SpikeTrains
    TT::Vector{Vector{Float64}}
    I::Vector{Int}
    rec_start::Float64  # usually just 0.0
    rec_finish::Float64 # defaults to maximum(map(maximum,TT))
    stim_start::Float64 # the start and end times of stimulus presentation
    stim_finish::Float64
    metadata::Dict{Any,Any}

    # A note on the field I : it is a list of labels. That is, the label for
    # TT[i] is I[i]. I will make sure this is enforced throughout the rest of
    # this project.

    function SpikeTrains(TT, I=1:length(TT);
        stim_start=0.0, stim_finish=maximum(map(maximum, TT)),
        rec_start=0.0, rec_finish=max(stim_finish, maximum(map(maximum, TT))),
        kwargs...)
        # use I to select which spike trains to keep, and while we're at it make
        # sure that TT ends up being the right type.
        trains = [map(Float64, TT[i][rec_start .< TT[i] .< rec_finish]) for i in I]
        new(trains, I, rec_start, rec_finish, stim_start, stim_finish, Dict(kwargs))
    end
end

n_cells(ST::SpikeTrains) = length(ST.I)

function show(io::IO, ST::SpikeTrains)
    println(io, "Spike times for $(n_cells(ST)) neurons.")
    println(io, "Recording start/finish: $(ST.rec_start) / $(ST.rec_finish) ")
    println(io, "Stimulus start/finish:  $(ST.stim_start) / $(ST.stim_finish)")
    show_metadata(io,ST)
end

"""
    histogram(ST, binsize; t0=ST.t_i, tmax=ST.t_f, N_bins=ceil(Int,(tmax - t0)/binsize))::Matrix{Float64}

Returns a `N_bins` by `N_cells` matrix where the `i,j`th entry is the number of
spikes that occurred in bin `i` from neuron `j`. By default returns a histogram
of spikes that happened during stimulus presentation. To truncate/extend the
histogram, specify `N_bins`. Truncation ignores spikes that happen later than
the end of the bins; extension pads with 0s.

"""
function histogram(ST::SpikeTrains, binsize::Float64;
    t0::Float64=ST.stim_start, tmax::Float64=ST.stim_finish,
    N_bins::Int=ceil(Int,(tmax - t0)/binsize))

    N_cells = n_cells(ST)
    # N_bins = ceil(Int,(tmax - t0)/binsize)
    spike_hist = zeros(N_bins, N_cells)
    for (j,S) = enumerate(ST.TT)
        Stilde = ceil.(Int, S[t0 .< S .< tmax]/binsize)
        max_idx = max(maximum(Stilde), N_bins)
        spike_hist[:,j] = sparsevec(Stilde, ones(length(Stilde)), max_idx)[1:N_bins]
        # sparsevec has a nice constructor that automatically does the counting
        # I want, so...
        # spike_hist[:,j] = sparsevec(ceil(Int, Stilde/binsize), ones(length(Stilde)), N_bins)
    end
    return spike_hist
end

"""
    raster(ST, binsize; t0, tmax)::BitMatrix

Returns `histogram(ST, binsize; t0=t0, tmax=tmax) .> 0`
"""
raster(ST::SpikeTrains, binsize::Float64; t0::Float64=ST.stim_start, tmax::Float64=ST.stim_finish, N_bins::Int=ceil(Int,(tmax - t0)/binsize)) = histogram(ST, binsize; t0=t0, tmax=tmax, N_bins=N_bins) .> 0
