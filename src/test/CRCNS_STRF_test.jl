using MAT
include("../probability/probability.jl")
include("../spikes/spikes.jl")
include("../stimulus/stimulus.jl")
using Probability, Spikes, Stimulus
include("../util/misc.jl")
include("../util/nonlinearities.jl")

data_dir = "Data"
fname = "20080516_R1.mat"
rec_idx = 1
println("* Loading a CRCNS matlab file: $fname")
println("* Creating CRCNS_Stimulus object.")
stim = CRCNS_Stimulus(joinpath(Stimulus.default_CRCNS_dir, data_dir, fname), rec_idx; verbose=true)
println(stim)
println("* Creating SpikeTrains object.")
spikes = Spikes.CRCNS_get_spikes_from_file(joinpath(Stimulus.default_CRCNS_dir, data_dir, fname), rec_idx)
println(spikes)

println("* Computing spike histograms with bin size = $(frame_time(stim))")
spike_hist = histogram(spikes, frame_time(stim))

println("* Computing STRFs")
STRFs = compute_STRFs(spike_hist, stim)

using PyPlot

println("* Creating plots...")
for cell = 1:7
    fig = figure("Cell $cell, STRF")
    # imshow(STRFs[:,cell,:], aspect="auto", cmap="gray")
    imshow(matrix_form(STRFs[cell]), aspect="auto", cmap="gray")
    colorbar()
    xlabel("time")
    ylabel("x coordinate")
end

# Now, let's simulate a neuron that is driven only by its receptive field
println("* Beginning simulation: convolving STRF 1 with stimulus")
r1, tau1 = STRF_response(STRFs[1], stim, flip_STRF_time=true)
n1 = spike_hist[:,1]

Q(u::Vector,v::Vector) = norm(u-v) / length(u)

c_range = (1.0:0.1:maximum(n1)) / tau1
h_range = [10.0^k for k in 2.0:0.1:4.0]
x0_range = decimal_round(minimum(r1),2):0.001:decimal_round(maximum(r1),2)
Q_vals = zeros(length(c_range), length(h_range), length(x0_range))
theta_ranges = [c_range, h_range, x0_range]
println("* Scaling response")
L1, theta1, Q1 = scale_response(r1, n1, sigmoid, Q; d=3, ranges=theta_ranges, verbose=true, save_fun=Q_vals)
