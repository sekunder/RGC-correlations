"""
    decimal_round(x, k)

    rounds `x` away from 0 to the `k`th place
"""
decimal_round(x::Float64, k::Integer) = x > 0 ? ceil(x * 10.0^k) / 10.0^k : floor(x * 10.0^k) / 10.0^k
