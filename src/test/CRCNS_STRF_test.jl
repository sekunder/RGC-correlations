
include("../probability/probability.jl")
include("../spikes/spikes.jl")
include("../stimulus/stimulus.jl")
using Probability, Spikes, Stimulus

data_dir = "Data"
fname = "20080516_R1.mat"
rec_idx = 1
println("Loading a CRCNS matlab file: $fname")
println("Creating CRCNS_Stimulus object.")
CS = CRCNS_Stimulus(joinpath(Stimulus.default_CRCNS_dir, data_dir, fname), rec_idx; verbose=true)
