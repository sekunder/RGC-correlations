
"""
    ReliableInformationDistribution

Distribution based on Schneidman et al 2011, where only "reliable" patterns are
used to influence the coefficients in the max-ent distribution.

"""
type ReliableInformationDistribution
    X::BitMatrix
    I::Vector{Int}
    metadata::Dict{Any,Any}
    cache::Dict{Any,Any}

    function ReliableInformationDistribution(X::Union{Matrix{Bool},BitMatrix}, I=1:size(X,1); kwargs)
        # MAYBEDO implement this if it seems necessary
    end
end
