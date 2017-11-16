
include("../util/init.jl")
# include("../util/metadata.jl")
include("../util/CRCNS/analysis_functions.jl")
using JLD, PyPlot

# plan: collect entropy values. 
D = CRCNS_collect_entropy(;verbose=2)
