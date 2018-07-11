
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
        hidden_keys = get(A.metadata, :hidden, Any[])
        for k in setdiff(keys(A.metadata), hidden_keys)
            println(io, "\t$k : $(A.metadata[k])")
        end
    end
end

"""
    hide_metadata!(A, k)

"Hides" metadata with key `k`. This means that `show_metadata` will only show
that this key exists, but not the value stored for it. This could be useful if,
say, you want to store the spike quality in a `SpikeTrains` object, which would
be a very long vector. If `k` is not a key in `A.metadata`, has no effect.

"""
function hide_metadata!(A, k)
    if haskey(A.metadata, k)
        A.metadata[:hidden] = union(get!(A.metadata,:hidden,Any[]), [k])
    end
end

"""
    unhide_metadata!(A, k)

"Unhides" metadata with key `k`. This means that `show_metadata` will show both
the key and the value for `k`. If `k` is not a key in `A.metadata`, has no
effect.

"""
function unhide_metadata!(A, k)
    if haskey(A.metadata, k)
        A.metadata[:hidden] = setdiff(get!(A.metadata,:hidden,Any[]), [k])
    end
end
