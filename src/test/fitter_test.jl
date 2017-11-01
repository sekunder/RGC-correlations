include("../probability/probability.jl")
include("../spikes/spikes.jl")
include("../stimulus/stimulus.jl")
using Probability, Spikes, Stimulus

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
spikes = Spikes.CRCNS_get_spikes_from_file(joinpath(Stimulus.default_CRCNS_dir, data_dir, fname), rec_idx)
println(spikes)

println("Loading stimulus cause of poor planning")
stim = CRCNS_Stimulus(joinpath(Stimulus.default_CRCNS_dir, data_dir, fname), rec_idx; verbose=true)
println("Computing spike histogram and raster with bin size = $(frame_time(stim))")
X = transpose(raster(spikes, frame_time(stim)))

println("Data distribution (i.e. histogram of patterns)")
tic()
P_data = data_model(X)
toc()
println(P_data)

println("First order model (i.e. Bernoulli code)")
tic()
P_1 = first_order_model(X)
toc()
println(P_1)

println("Second order model (i.e. Ising model)")
tic()
# P_2 = second_order_model(X, verbose=true, algorithm=:LD_MMA, maxeval=20) #algorithm=:LD_MMA,
P_2 = second_order_model(X; verbose=true, more_verbose=true)
# (F_opt, J_opt, stop, Jseed, mu, F_used) = second_order_model(X; verbose=true, more_verbose=true, print_eval=10, maxeval=100)
toc()
println(P_2)

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
