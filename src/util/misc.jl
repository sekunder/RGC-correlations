"""
    decimal_round(x, k)

    rounds `x` away from 0 to the `k`th place
"""
decimal_round(x::Float64, k::Integer) = x > 0 ? ceil(x * 10.0^k) / 10.0^k : floor(x * 10.0^k) / 10.0^k

"""
    remove_extension(fname)

Removes the trailing ".ext" from `fname`
"""
remove_extension(fname) = join(split(fname,'.')[1:end-1],'.')
