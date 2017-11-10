
"""
    show_metadata(io, A::Any)

Displays the entries in the `metadata` field of `A`. Since nearly every type
I've written for this project involves a metadata field, I figured it would be
beneficial to move this to a universally accessible location.

"""
function show_metadata(io::IO, A::Any)
    if isempty(A.metadata)
        println(io, "No metadata found")
    else
        println(io, "Metadata:")
        for (k,v) in A.metadata
            println(io, "\t$k : $v")
        end
    end
end

#TODO(big) either scrap existing data, or write some compatibility scripts to
#handle the fact that this could break some JLD files (not that I've actually
#made any yet, but it's the kind of thing to keep in mind...)

#TODO(small) add "hidden metadata", i.e. information where I only want to know that a key exists
#TODO(small) add to Bernoulli
#TODO(small) add to Data
#TODO(small) add to Ising
#TODO(small) add to Spikes
#TODO(small) add to Stimulus
#TODO(small) hide_metadata! and unhide_metadata!

"""
    metadata(A::Any, k::Any)

Convenience function, returns `A`'s metadata value for key `k`, or `:none` if it
is not set.
"""
metadata(A::Any, k::Any) = get(A.metadata, k, :none)

"""
    metadata!(A::Any, k::Any, v::Any=:none)

Convenience function, sets metadata for `A`. Can be called without third
argument to ensure key `k` exists in `A`'s metadata.
"""
metadata!(A::Any, k::Any, v::Any=:none) = get!(A.metadata, k, v)
