include("../probability/probability.jl")
include("../spikes/spikes.jl")
include("../stimulus/stimulus.jl")
using Probability, Spikes, Stimulus
include("../util/constants.jl")

# 1. load up the CRCNS spikes data
# 2. get the spike raster
# 3. create data model, first model, second model
# 4. compute expectation matrices of each
# 5. plot comparisons

data_dir = "Data"
fname = "20080516_R1.mat"
rec_idx = 1
println("Loading a CRCNS matlab file: $fname")
println("Creating SpikeTrains object.")
spikes = Spikes.CRCNS_get_spikes_from_file(joinpath(default_CRCNS_dir, data_dir, fname), rec_idx)
println(spikes)

println("Loading stimulus cause of poor planning")
stim = CRCNS_Stimulus(joinpath(default_CRCNS_dir, data_dir, fname), rec_idx; verbose=true)
println("Computing spike histogram and raster with bin size = $(frame_time(stim))")
X = transpose(raster(spikes, frame_time(stim)))
N_neurons,N_samples = size(X)

println("--------------------------------------------------------------------------------")
println("* Second order model (i.e. Ising model). Two different fitters: gradient_optimizer, NLopt using LD_LBFGS")
println("* Also two different functions: L_X and K_X")
println("--------------------------------------------------------------------------------")
println("* Varying the learning rate parameter for gradient_optimizer")
println("\tThe goal is to decide what the default algorithm should be.")
lr_options = [1.0, 1.0 * N_neurons, 10.0, 1.0*N_neurons^2]
L_naive = Dict()
for lr in lr_options
    tic()
    P_2 = second_order_model(X; verbose=2, algorithm=:naive, lr=lr)
    t = toc()
    L_naive[lr] = (P_2, t, P_2.metadata[:opt_val])
end
println("--------------------------------------------------------------------------------")
println("* Now using NLopt with LD_LBFGS algorithm")
tic()
P_2_NLopt_L = second_order_model(X; verbose=true, algorithm=:LD_LBFGS)
t_NLopt_L = toc()
L_NLopt = (P_2_NLopt_L, t_NLopt_L, P_2_NLopt_L.metadata[:opt_val])
println("--------------------------------------------------------------------------------")

println("* Using gradient_optimizer, here are the function values and times (followed by learning rate)")
for (k,v) in L_naive
    println("$(v[2])\t$(v[3]) s\t(lr = $k)")
end
println("* Compare to NLopt:")
println("$(L_NLopt[2])\t$(L_NLopt[3]) s")

# K_naive = Dict()
# for lr in lr_options
#     tic()
#     P_2 = second_order_model(X; verbose=2, algorithm=:naive, lr=lr, force_MPF=true, print_eval=1, maxeval=10)
#     t = toc()
#     K_naive[lr] = (P_2,t, P_2.metadata[:opt_val])
# end
#
# tic()
# P_2_NLopt_K = second_order_model(X; verbose=true, algorithm=:LD_LBFGS, force_MPF=true)
# t_NLopt_K = toc()
# K_NLopt = (P_2_NLopt_K, t_NLopt_K, P_2_NLopt_K.metadata[:opt_val])


# # MANUAL GRADIENT DESCENT
# N_neurons, N_samples = size(X)
# F_X(J,g) = loglikelihood(X, J, g; mu_X =X * X' / N_samples)
#
# begin
#     n_attempts=100
#     print_attempts=10
#     lr = 10.0; whoopscount=0;
#     println("Gonna just do gradient descent for $n_attempts steps with a learning rate")
#     Jseed = rand(N_neurons,N_neurons); Jseed = (Jseed + Jseed') / (2 * N_neurons)
#     J_0 = Jseed[:]
#     g = zeros(N_neurons^2)
#
#     F_vals = zeros(0)
#     J_vals=Vector{Vector{Float64}}(0)
#     g_vals = Vector{Vector{Float64}}(0)
#
#     push!(J_vals, J_0)
#     push!(F_vals, F_X(J_0, g))
#     push!(g_vals, g)
#     max_F = F_vals[1]; max_F_idx = 1;
#     attempt = 1
#     println("Step\tF_X(J,θ)\tdF\t|gradient|\tlr")
#     println("$attempt\t$(F_vals[attempt])\t$(F_vals[attempt])\t$(sqrt(dot(g,g)))\t$lr")
#     attempt += 1
#     tic()
#     while attempt <= n_attempts
#         push!(J_vals, J_vals[attempt-1] + lr * g_vals[attempt-1])
#         push!(F_vals, F_X(J_vals[attempt], g))
#         push!(g_vals, g)
#         if mod(attempt, print_attempts) == 0
#             println("$attempt\t$(F_vals[attempt])\t$(F_vals[attempt]-F_vals[attempt-1])\t$(lr * sqrt(dot(g,g)))\t$lr")
#         end
#         if F_vals[attempt] < F_vals[attempt-1]
#             whoopscount += 1
#             # lr = 2.0 - exp(-whoopscount)
#             lr /= 2.0
#         end
#         attempt += 1
#     end
#     toc()
# end
