
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
    "GrayScaleStimuli" => "https://github.com/sekunder/GrayScaleStimuli.jl",
    # "PrettyLogging" => "https://github.com/sekunder/PrettyLogging.jl"
    )
for (pkg_name, pkg_url) in unregistered_packages
    if !(pkg_name in keys(Pkg.installed()))
        Pkg.clone(pkg_url)
    end
    # Pkg.update(pkg_name)
end
@everywhere using BinaryVectorProbability, Spikes, GrayScaleStimuli, PrettyLogging

@everywhere include("../util/constants.jl")
@everywhere include("../util/misc.jl")
@everywhere include("../util/metadata.jl")
@everywhere include("../util/process_args.jl")
@everywhere include("../util/database.jl")
