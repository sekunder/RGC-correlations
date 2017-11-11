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

"""
    random_subset(A,k)

Returns a random subset of k elements from collection A. Uses Reservoir
Sampling, see https://en.wikipedia.org/wiki/Reservoir_sampling#Algorithm_R

"""
function random_subset(A,k::Integer)
    n = length(A)
    # assert(n >= k, "Reservoir Sample error: To few items to sample from")
    R = Array(typeof(A[1]), k)
    R[1:k] = A[1:k]
    for i = (k+1):n
        j = rand(1:i)
        if j <= k
            R[j] = A[i]
        end
    end
return R
end

"""
    random_subset!(A,R)

Returns a random subset of length(R) elements from collection A, modifying R in
place. See `random_subset` for more details.

"""
function random_subset!(A, R::Vector)
    n = length(A)
    k = length(R)
    # assert(n >= k, "Reservoir Sample error: Reservoir too large")
    R[1:k] = A[1:k]
    for i = (k+1):n
        j = rand(1:i)
        if j <= k
            R[j] = A[i]
        end
    end
    return R
end
