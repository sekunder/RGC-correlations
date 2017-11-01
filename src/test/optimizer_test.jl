include("../util/optimizers.jl")

function G(x, g)
    if length(g) > 0
        g[1] = 2.0*x[1] - 1.0
    end
    return (x[1] - 3.0)*(x[1] + 2.0)
end

x0 = [5.0]
println("Attempting to optimize G(x) = x^2 - x - 6, starting from $x0")
println("First things first: G(5) = $(G(x0,[]))")
(F_opt, x_opt, stop) = gradient_optimizer(G, x0; verbose=2, print_eval=1, lr=0.2)

println("return valuse from grad_opt: $F_opt, $x_opt, $stop")
println("evaluating function at x_opt: $(G(x_opt,[]))")
