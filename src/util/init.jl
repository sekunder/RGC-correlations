
# These inclusions/etc. go at the start of every script I write, basically, so
# I'm just putting them in one place

include("../probability/probability.jl")
include("../spikes/spikes.jl")
include("../stimulus/stimulus.jl")
@everywhere using Spikes, Stimulus, Probability
include("../util/constants.jl")
include("../util/misc.jl")
include("../util/metadata.jl")
include("../util/process_args.jl")
