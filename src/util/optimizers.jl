
"""
    gradient_optimizer(obj_fun, x0; kwargs...)

Performs a very rudimentary version of gradient descent/ascent. Follows basic
pattern of NLopt package: `obj_fun` is a function which accepts arguments `x`, `g`; it
returns the objective function at `x` and sets the gradient `g` in place. The
value `x0` is an initial guess. Default behavior is to descend along the
gradient to minimize the function, setting keyword argument `objective` to
`:max` changes this behavior (default is `:min`). Returns a tuple
`(F_opt, x_opt, stop)` where `stop` is an array of the stopping criteria met.

List of keywords and default values:
 * `objective = :min` minimize `obj_fun`, setting to `:max` will maximize
 * `verbose = 0` Options are 0 (no output), 1 (some output), 2 (all the output)
 * `print_eval = 100` In maximum verbosity, waits this many evaluations before printing
 * `lr = 1.0` A "learning rate", the gradient is scaled by this much.
 * `adjust_lr = lr > 1` If true, every time a smaller value of `obj_fun` is reached, `lr` is halved.

Stopping criteria (except vanishing gradient) can be set with keyword arguments.
 * `maxeval = 1000` stop after evaluating `obj_fun` this many times
 * `xtol_abs = 0.0` stop when `|x_cur - x_prev| < xtol_abs`
 * `xtol_rel = 0.0` stop when `|x_cur - x_prev|/|x_cur| < xtol_rel`
 * `ftol_abs = 0.0` stop when `|F_cur - F_prev| < ftol_abs`
 * `ftol_rel = 0.0` stop when `|F_cur - F_prev|/|F_cur| < ftol_rel`

"""
function gradient_optimizer(obj_fun, x0; objective=:min,
    verbose::Int=0, print_eval::Int=100,
    lr::Float64=1.0, adjust_lr::Bool=(lr > 1.0),
    maxeval::Int=1000,
    xtol_abs::Float64=0.0, xtol_rel::Float64=0.0,
    ftol_abs::Float64=0.0, ftol_rel::Float64=0.0,
    kwargs...)

    # println("The type of obj_fun is $(typeof(obj_fun))")
    # obj_fun = (x,g) -> obj_fun(x,g)
    # println("Now the type of obj_fun is $(typeof(obj_fun))")
    # bookkeeping set up
    n_vars = length(x0)
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
    g_cur = zeros(x0); g_prev = zeros(x0); g_opt = zeros(x0)
    x_cur = zeros(x0); x_prev = zeros(x0); x_opt = zeros(x0)
    F_cur = 0.0; F_prev = 0.0; F_opt = -obj_sign * Inf;
    dF_cur = 0.0; dx_cur = 0.0; dx_next = 0.0
    eval = 0; eval_opt = 0
    x_prev[:] = x0[:]

    whoopscount = 0

    if verbose > 0
        println("gradient_optimizer: $(obj_sign > 0 ? "maximizing" : "minimizing") function of $n_vars variable$(n_vars > 1 ? "s" : "")")
    end
    if verbose > 1
        println("\tlr           : $lr ($(adjust_lr ? "will adjust" : "no adjust"))")
        println("\tmaxeval      : $maxeval")
        println("\txtol abs/rel : $xtol_abs / $xtol_rel")
        println("\tftol abs/rel : $ftol_abs / $ftol_rel")
        println()
        println("\teval\tF(x_cur)\t\tdF\t|dx_next|\tlr")
    end
    while stop.chunks[1] == 0
        # get current value of x_cur from x_prev and gradient
        x_cur[:] = x_prev[:] + obj_sign * lr * g_prev[:]
        # evaluate at new x_cur
        F_cur = obj_fun(x_cur, g_cur)
        dF_cur = F_cur - F_prev
        dx_next = lr * norm(g_cur) # size of next step

        # perform updates as necessary
        if F_cur * obj_sign > F_opt * obj_sign
            # if we're more extreme than current best, update
            g_opt[:] = g_cur[:]
            x_opt[:] = x_cur[:]
            F_opt = F_cur
            eval_opt = eval
        end
        if adjust_lr && (dF_cur * obj_sign < 0)
            whoopscount += 1
            lr /= 2
        end

        # print if necessary
        if verbose > 1 && mod(eval, print_eval) == 0
            println("\t$eval\t$F_cur\t$dF_cur\t$dx_next\t$lr$(eval_opt == eval ? " *" : "")")
        end

        # increment counter
        eval += 1

        #push cur into prev slot
        x_prev[:] = x_cur[:]
        g_prev[:] = g_cur[:]
        F_prev = F_cur


        # Check stopping conditions
        # 1 maxeval
        # 2 opt reached (dF = 0)
        # 3 xtol_abs reached (lr * norm(g_cur) < xtol_abs)
        # 4 xtol_rel reached (lr * norm(g_cur) / norm(x_cur) < xtol_rel)
        # 5 ftol_abs reached (dF_cur < ftol_abs)
        # 6 ftol_rel reached (|dF_cur / F_cur| < ftol_rel )
        Fstep_rel = F_cur == 0 ? Inf : abs(dF_cur / F_cur)
        xstep_rel = norm(x_cur) == 0 ? Inf : abs(dx_next / norm(x_cur))
        stop = stop | [eval >= maxeval,
                        dF_cur == 0.0,
                        dx_next < xtol_abs,
                        xstep_rel < xtol_rel,
                        abs(dF_cur) < ftol_abs,
                        Fstep_rel < ftol_rel]
    end
    if verbose > 0
        println("gradient_optimizer: Stopping criteria reached: $(stopping_criteria[stop])")
        println("gradient_optimizer: optimal value found: $F_opt")
    end
    if verbose > 1
        println("gradient_optimizer: opt. value reached at evaluation $eval_opt / $eval")
        println("gradient_optimizer: Learning rate adjusted $whoopscount time$(whoopscount != 1 ? "s" : "")")
        println("gradient_optimizer: final learning rate: $lr")
    end

    return F_opt, x_opt, stop
end
