include("../probability/probability.jl")
include("../spikes/spikes.jl")
include("../stimulus/stimulus.jl")
using Probability, Spikes, Stimulus
using PyPlot
include("../util/constants.jl")

# 1. load up the CRCNS spikes data
# 2. get the spike raster
# 3. create data model, first model, second model
# 4. compute expectation matrices of each
# 5. plot comparisons

################################################################################
#### SET UP
################################################################################
begin
    data_dir = "Data"
    fname = "20080516_R1.mat"
    rec_idx = 1
    println("* Loading a CRCNS matlab file: $fname")
    println("* Creating SpikeTrains object.")
    spikes = Spikes.CRCNS_get_spikes_from_file(joinpath(CRCNS_dir, data_dir, fname), rec_idx)
    println(spikes)

    println("* Loading stimulus cause of poor planning")
    stim = CRCNS_Stimulus(joinpath(CRCNS_dir, data_dir, fname), rec_idx; verbose=true)
    println("* Computing spike histogram and raster with bin size = $(frame_time(stim))")
    X = transpose(raster(spikes, frame_time(stim)))
    N_neurons,N_samples = size(X)
end

################################################################################
#### OBJECTIVE FUNCTIONS
################################################################################
mu_data = X * X' / size(X,2)
L_X(J,g) = loglikelihood(X, J, g; mu_X=mu_data)
K_X(J,g) = MPF_objective(X, J, g)

################################################################################
#### FIT MODELS
################################################################################

begin # Data model
    println("* Fitting data model (i.e. histogram of patterns)")
    tic()
    P_X = data_model(X)
    toc()
    println(P_X)
end

begin # first order model
    println("* Fitting first order model (i.e. Bernoulli Code)")
    tic()
    P_1 = first_order_model(X)
    toc()
    println(P_1)
end

begin # second order, L
    println("* Fitting second order model (i.e. Ising model), using loglikelihood")
    tic()
    P_2_L = second_order_model(X; verbose=true)
    toc()
    println(P_2_L)
end

# begin # second order, K
#     println("* Fitting second order model using MPF")
#     tic()
#     P_2_K = second_order_model(X; verbose=true, force_MPF=true)
#     toc()
#     println(P_2_K)
# end

################################################################################
#### EXPECTATION MATRICES
################################################################################
mu_X = expectation_matrix(P_X)
mu_1 = expectation_matrix(P_1)
mu_2_L = expectation_matrix(P_2_L)
# mu_2_K = expectation_matrix(P_2_K)

################################################################################
#### PLOTS
################################################################################

begin # Scatter P_* vs. P_X using log-log plot
    scatter_fig = figure("P_1, P_2 vs. P_X; $fname / 1", tight_layout=true)

    plot([0,1],[0,1])
    scatter(get_pdf(P_X), get_pdf(P_1), label=L"P_1", color="blue", s=5)
    scatter(get_pdf(P_X), get_pdf(P_2_L), label=L"P_2(L)", color="red", s=5)
    # scatter(get_pdf(P_X), get_pdf(P_2_K), label=L"P_2(K)", color="green", s=5)
    title("First Order Model vs. Data")
    xlabel(L"P_X")
    ylabel(L"P_1")
    # xscale("log")
    # yscale("log")
    # xlim(1e-12,maximum(get_pdf(P_X)))
    # ylim(1e-15,maximum(get_pdf(P_1)))
    xlim(0,maximum(get_pdf(P_X)))
    ylim(0,max(maximum(get_pdf(P_1)), maximum(get_pdf(P_2_L))))
    legend(loc="upper left")
end

################################################################################
#### SAMPLING
################################################################################
samples = Dict()
begin # sample from P_X
    println("* Sampling from data distribution P_X")
    tic()
    X_X = random(P_X, N_samples)
    toc()
    samples["P_X"] = X_X
end

begin # sample from P_1
    println("* Sampling from Bernoulli distribution P_1")
    tic()
    X_1 = random(P_1, N_samples)
    toc()
    samples["P_1"] = X_1
end

begin # sample from P_2_L, using exact
    println("* Sampling from P_2_L, using exact sampling")
    tic()
    X_2_L_exact = random(P_2_L, N_samples)
    toc()
    samples["P_2_L_exact"] = X_2_L_exact
end

begin #sample from P_2_L, using gibbs
    println("* Sampling from P_2_L, using gibbs sampling")
    tic()
    X_2_L_gibbs = random(P_2_L, N_samples, true)
    toc()
    samples["P_2_L_gibbs"] = X_2_L_gibbs
end

################################################################################
#### MORE PLOTS
################################################################################
for (name, sample) in samples
    mu_sample = sample * sample'/N_samples
    figure("$name compared to real data", tight_layout=true)
    subplot(131)
    imshow(mu_sample, interpolation="none")
    colorbar()
    title("$name")

    subplot(132)
    imshow(mu_data, interpolation="none")
    colorbar()
    title("Real data")

    subplot(133)
    imshow(mu_sample-mu_data, interpolation="none")
    colorbar()
    title("difference")
end
