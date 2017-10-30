using MAT
include("../probability/probability.jl")
include("../spikes/spikes.jl")
include("../stimulus/stimulus.jl")
using Probability, Spikes, Stimulus

data_dir = "Data"
fname = "20080516_R1.mat"
rec_idx = 1
println("Loading a CRCNS matlab file: $fname")
println("Creating CRCNS_Stimulus object.")
stim = CRCNS_Stimulus(joinpath(Stimulus.default_CRCNS_dir, data_dir, fname), rec_idx; verbose=true)
println(stim)
println("Creating SpikeTrains object.")
spikes = Spikes.CRCNS_get_spikes_from_file(joinpath(Stimulus.default_CRCNS_dir, data_dir, fname), rec_idx)
println(spikes)

println("Computing spike histograms with bin size = $(frame_time(stim))")
spike_hist = histogram(spikes, frame_time(stim))

println("Computing STRFs")
STRFs = compute_STRFs(spike_hist, stim)

using PyPlot

for cell = 1:7
    fig = figure("Cell $cell, STRF")
    # imshow(STRFs[:,cell,:], aspect="auto", cmap="gray")
    imshow(matrix_form(STRFs[cell]), aspect="auto", cmap="gray")
    colorbar()
    xlabel("time")
    ylabel("x coordinate")
end
