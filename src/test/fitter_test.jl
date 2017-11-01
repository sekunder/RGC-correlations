include("../probability/probability.jl")
include("../spikes/spikes.jl")
include("../stimulus/stimulus.jl")
using Probability, Spikes, Stimulus
using PyPlot

# 1. load up the CRCNS spikes data
# 2. get the spike raster
# 3. create data model, first model, second model
# 4. compute expectation matrices of each
# 5. plot comparisons

################################################################################
#### SET UP
################################################################################
data_dir = "Data"
fname = "20080516_R1.mat"
rec_idx = 1
println("* Loading a CRCNS matlab file: $fname")
println("* Creating SpikeTrains object.")
spikes = Spikes.CRCNS_get_spikes_from_file(joinpath(Stimulus.default_CRCNS_dir, data_dir, fname), rec_idx)
println(spikes)

println("* Loading stimulus cause of poor planning")
stim = CRCNS_Stimulus(joinpath(Stimulus.default_CRCNS_dir, data_dir, fname), rec_idx; verbose=true)
println("* Computing spike histogram and raster with bin size = $(frame_time(stim))")
X = transpose(raster(spikes, frame_time(stim)))

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
#     println(P_2_L)
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
    scatter(get_pdf(P_X), get_pdf(P_1), label="P_1", color="blue", s=5)
    scatter(get_pdf(P_X), get_pdf(P_2_L), label="P_2", color="red", s=5)
    title("First Order Model vs. Data")
    xlabel(L"P_X")
    ylabel(L"P_1")
    # xscale("log")
    # yscale("log")
    xlim(1e-12,maximum(get_pdf(P_X)))
    ylim(1e-15,maximum(get_pdf(P_1)))
    legend(loc="upper left")
end



################################################################################
#### CLEAN UP
################################################################################
close_all_the_figs = false
if close_all_the_figs
    map(close,[scatter_fig])
end
