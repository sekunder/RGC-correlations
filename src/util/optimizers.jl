
"""
    gradient_optimizer(F, x0; kwargs...)

Performs a very rudimentary version of gradient descent/ascent. Follows basic
pattern of NLopt package: `F` is a function which accepts arguments `x`, `g`; it
returns the objective function at `x` and sets the gradient `g` in place. The
value `x0` is an initial guess. Default behavior is to descend along the
gradient to minimize the function, setting keyword argument `objective` to
`:max` changes this behavior (default is `:min`).

List of keywords and default values:
 * `objective = :min` minimize `F`, setting to `:max` will maximize
 * `verbose = 0` Options are 0 (no output), 1 (some output), 2 (all the output)
 * `print_eval = 100` In maximum verbosity, waits this many evaluations before printing
 * `lr = 1.0` A "learning rate", the gradient is scaled by this much.
 * `adjust_lr = false` If true, every time a smaller value of `F` is reached, `lr` is halved.

Stopping criteria are also passed as keyword arguments:
 * `maxeval = 1000` maximum number of function evaluations
 * `xtol_abs = 0.0` stop when `|x_cur - x_prev| < xtol_abs`
 * `xtol_rel = 0.0` stop when `|x_cur - x_prev|/|x_cur| < xtol_rel`
 * `ftol_abs = 0.0` stop when `|F_cur - F_prev| < ftol_abs`
 * `ftol_rel = 0.0` stop when `|F_cur - F_prev|/|F_cur| < ftol_rel`

"""
function gradient_optimizer(F, x0; objective=:min,
    verbose=0, print_eval=100, lr=1.0, adjust_lr=(lr > 0.0),
    maxeval=1000, xtol_abs = 0.0, xtol_rel = 0.0, ftol_abs = 0.0, ftol_rel = 0.0,
    kwargs...)

    # bookkeeping set up
    dkwargs = Dict(kwargs)
    obj_sign = (objective == :min ? -1 : 1)
    stopping_criteria = [:maxeval_reached,
                        :opt_reached,
                        :xtol_abs_reached,
                        :xtol_rel_reached,
                        :ftol_abs_reached,
                        :ftol_rel_reached]
    stop = falses(stopping_criteria)

    # gradient algorithm set up
    g = zeros(x0); g_prev = zeros(x0)
    x_prev = zeros(x0); x_prev[:] = x0[:]
    x = zeros(x0); dx_cur = 0

    F_cur = 0.0; F_prev = 0.0
    dF_cur = 0.0
    F_opt = -obj_sign * Inf; eval_opt = 0;
    x_opt[:] = zeros(x0); g_opt = zeros(x0)
    eval = 0
    whoopscount = 0

    if verbose > 1
        println("Eval\tF(x)\t\tdF\t|gradient|\tlr")
    end
    while stop.chunks[1] == 0
        # get current value of x from previous x and gradient
        x[:] = x_prev[:] + obj_sign * lr * g_prev[:]
        # evaluate at new x
        F_cur = F(x, g)
        dF_cur = F_cur - F_prev
        dx_next = lr * norm(g) # size of next step
        if verbose > 1 && mod(eval, print_eval) == 0
            println("$eval\t$F_cur\t$dF_cur\t$(lr * sqrt(dot(g,g)))\t$lr")
        end

        # perform updates as necessary
        # if (F_cur - F_opt) * obj_sign > 0
        if F_cur * obj_sign < F_opt * obj_sign
            # if we're more extreme than current best, update
            F_opt = F_cur
            eval_opt = eval
        end
        if adjust_lr && (dF_cur * obj_sign < 0)
            whoopscount += 1
            lr /= 2
        end
        eval += 1

        #push cur into prev slot
        x_prev[:] = x[:]
        g_prev[:] = g[:]
        F_prev = F_cur


        # Check stopping conditions
        # 1 maxeval
        # 2 opt reached (dF = 0)
        # 3 xtol_abs reached (lr * norm(g) < xtol_abs)
        # 4 xtol_rel reached (lr * norm(g) / norm(x) < xtol_rel)
        # 5 ftol_abs reached (dF_cur < ftol_abs)
        # 6 ftol_rel reached (|dF_cur / F_cur| < ftol_rel )
        Fstep_rel = F_cur == 0 ? Inf : abs(dF_cur / F_cur)
        xstep_rel = norm(x) == 0 ? Inf : abs(dx_next / norm(x))
        stop = stop | [eval >= maxeval,
                        dF_cur == 0.0,
                        dx_next < xtol_abs,
                        xtep_rel < xtol_rel,
                        dF_cur < ftol_abs,
                        Fstep_rel < ftol_rel]
    end
    if verbose > 0
        # TODO print out information about the optimal value, that sort of thing
        println("gradient_optimizer: Stopping criteria reached: $(stopping_criteria[stop])")
        println("gradient_optimizer: optimal value found: $F_opt")
        println("gradient_optimizer: opt. value reached at evaluation $eval_opt")
    end
    if verbose > 1
        println("gradient_optimizer: Learning rate adjusted $whoopscount time$(whoopscount > 1 ? "s" : "")")
        println("gradient_optimizer: final learning rate: $lr")
    end
end
