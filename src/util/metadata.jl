
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
