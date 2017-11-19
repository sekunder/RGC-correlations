
"""
    process_args(args; kwargs...)

Convenience function. Returns a `Dict{String,Any}` where `D[k]` is the value of
flag `-k` passed at the command line:
`--flag v` is returned as `D[flag] = v`
`-flag` is returned as `D[flag] = true`
The special key "0" returns all non-flag values in the order they appeared.
Note that "v" is always a bool flag and "verbose" is a parse flag with default
value 0

Keyword arguments:
 * `parse_flags::Vector{String}` Looks for commandline flags of the form "--flag v" and attempts to parse the value of `v` using `parse`
 * `parse_defaults::Dict{String,Any}` Same as above, but will return the default value specified if the flag is not found
 * `bool_flags::Vector{String}` Any flag of the form `-flag` will get picked up as `D["flag"] = true`, but adding it to `bool_flags` will return false if it's not there

"""
function process_args(args; bool_flags=String[], parse_flags=String[],
    parse_defaults=Dict{String,Any}(), kwargs...)

    D = Dict{String,Any}()
    dkwargs = Dict(kwargs)
    # bool_flags = ["v"; pop!(dkwargs, :bool_flags, String[])
    # parse_defaults = pop!(dkwargs, :parse_defaults, Dict{String,Any}())
    # # parse_flags = union(pop!(dkwargs, :parse_flags, ["n_trials", "verbose"]), keys(parse_defaults))
    # parse_flags = pop!(dkwargs, :parse_flags, ["verbose"])
    push!(bool_flags, "v")
    # push!(parse_flags, "")
    get!(parse_defaults, "verbose", 0)
    # for flag in parse_flags
    #     flag_indexes = find(x -> x == "--$flag", args)
    #     if !isempty(flag_indexes)
    #         D[flag] = parse(args[flag_indexes[end]+1])
    #         drop_indexes = sort([flag_indexes; flag_indexes + 1])
    #         args = args[setdiff(linearindices(args), drop_indexes)]
    #     elseif haskey(parse_defaults, flag)
    #         D[flag] = parse_defaults[flag]
    #     end
    # end
    # for flag in bool_flags
    #     flag_indexes = find(x -> x == "-$flag", args)
    #     D[flag] = !isempty(flag_indexes)
    #     args = args[setdiff(linearindices(args), flag_indexes)]
    # end
    idx = 1
    D["0"] = String[]
    while idx <= length(args)
        if startswith(args[idx], "--")
            D[args[idx][3:end]] = idx == length(args) ? true : args[idx+1]
            idx += 2
        elseif startswith(args[idx], "-")
            D[args[idx][2:end]] = true
            idx += 1
        else
            push!(D["0"], args[idx])
            idx += 1
        end
    end
    for flag in bool_flags
        get!(D, flag, false)
    end
    for (flag, val) in parse_defaults
        if haskey(D, flag)
            D[flag] = parse(D[flag])
        else
            D[flag] = val
        end
    end
    for flag in parse_flags
        if haskey(D, flag)
            D[flag] = parse(D[flag])
        end
    end
    return D
end
