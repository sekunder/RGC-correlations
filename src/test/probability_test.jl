include("../probability/probability.jl")
using Probability
N_small = 10; N_large =25; N_samples = 10000;

println("Creating binary white noise data; $N_large neurons, $N_samples samples")

X = rand(N_large, N_samples) .< 0.5
Ps = rand(N_large)
J_small = rand(N_small, N_small) / N_small^2
J_large = rand(N_large, N_large) / N_large^2

DD_small = DataDistribution(X, 1:10)
println("Created data distribution from first 10 neurons")
println(DD_small)


ID_small = IsingDistribution(J_small; comment="small J")
ID_large = IsingDistribution(J_large; comment="large J")

X_Ising_small = random(ID_small, 100)
