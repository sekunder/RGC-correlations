include("../probability/probability.jl")
include("../spikes/spikes.jl")
include("../stimulus/stimulus.jl")
using Probability, Spikes, Stimulus

# 1. load up the CRCNS spikes data
# 2. get the spike raster
# 3. create data model, first model, second model
# 4. compute expectation matrices of each
# 5. plot comparisons

data_dir = "Data"
fname = "20080516_R1.mat"
rec_idx = 1
println("Loading a CRCNS matlab file: $fname")
println("Creating SpikeTrains object.")
spikes = Spikes.CRCNS_get_spikes_from_file(joinpath(Stimulus.default_CRCNS_dir, data_dir, fname), rec_idx)
println(spikes)

println("Loading stimulus cause of poor planning")
stim = CRCNS_Stimulus(joinpath(Stimulus.default_CRCNS_dir, data_dir, fname), rec_idx; verbose=true)
println("Computing spike histogram and raster with bin size = $(frame_time(stim))")
X = transpose(raster(spikes, frame_time(stim)))

println("Data distribution (i.e. histogram of patterns)")
tic()
P_data = data_model(X)
toc()
println(P_data)

println("First order model (i.e. Bernoulli code)")
tic()
P_1 = first_order_model(X)
toc()
println(P_1)

println("Second order model (i.e. Ising model)")
tic()
P_2 = second_order_model(X)
toc()
println(P_2)
