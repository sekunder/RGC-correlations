
# These inclusions/etc. go at the start of every script I write, basically, so
# I'm just putting them in one place

# include("../probability/probability.jl")
# include("../spikes/spikes.jl")
# include("../stimulus/stimulus.jl")
# using Spikes, Stimulus, Probability

# Each of the above files has been moved into a separate git repo that the Julia
# package manager can use. Should allow for better parallel processing,
# hopefully?
unregistered_packages = Dict(
    "BinaryVectorProbability" => "https://github.com/sekunder/BinaryVectorProbability.jl.git",
    "Spikes" => "https://github.com/sekunder/Spikes.jl",
    "GrayScaleStimuli" => "https://github.com/sekunder/GrayScaleStimuli.jl"
)
for pkg_name, pkg_url in unregistered_packages
    if !(pkg_name in keys(Pkg.installed()))
        Pkg.clone(pkg_url)
    end
end
using BinaryVectorProbability, Spikes, GrayScaleStimuli

include("../util/constants.jl")
include("../util/misc.jl")
include("../util/metadata.jl")
include("../util/process_args.jl")
