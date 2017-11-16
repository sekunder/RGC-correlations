
"""
    process_args(args)

Convenience function. Returns a `Dict{String,Any}` where `D[k]` is the value of
flag `-k` passed at the command line (i.e., in `args`, one entry is `-k` and the
immediate next entry is stored in `D[k]`). Does some parsing for types and
allows for boolean flags. The special key "0" returns all non-flag values in the
order they appeared.

"""
function process_args(args)
    D = Dict{String,Any}()
    parse_flags = ["n_trials", "verbose"]
    bool_flags = ["v"]
    for flag in parse_flags
        flag_indexes = find(x -> x == "--$flag", args)
        if !isempty(flag_indexes)
            D[flag] = parse(args[flag_indexes[end]+1])
            drop_indexes = sort([flag_indexes; flag_indexes + 1])
            args = args[setdiff(linearindices(args), drop_indexes)]
        end
    end
    for flag in bool_flags
        flag_indexes = find(x -> x == "-$flag", args)
        D[flag] = !isempty(flag_indexes)
        args = args[setdiff(linearindices(args), flag_indexes)]
    end
    idx = 1
    D["0"] = String[]
    while idx <= length(args)
        if startswith(args[idx], "--")
            D[args[idx][3:end]] = args[idx+1]
            idx += 2
        elseif startswith(args[idx], "-")
            D[args[idx][2:end]] = true
            idx += 1
        else
            push!(D["0"], args[idx])
            idx += 1
        end
    end
    return D
end
