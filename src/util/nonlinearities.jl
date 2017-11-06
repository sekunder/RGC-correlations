# Some paramaterized nonlinearities with typical default values
"""
    sigmoid(x, theta)

Parameterized sigmoid.
`f(x; θ) = θ_1 / (1 + exp(-θ_2 * (x - θ_3)))`
"""
sigmoid(x, theta::Vector{Float64}=[1.0,1.0,0.0]) = theta[1] ./ (1 + exp(-theta[2] .* (x .- theta[3])))
