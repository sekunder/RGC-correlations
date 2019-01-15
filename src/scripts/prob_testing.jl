include("../util/init.jl")

using CSV

df = load_prob_db(CRCNS_db_prob)

df[:,:I_2] = (df[:,:H_1] .- df[:,:H_2]) ./ (df[:,:H_1] .- df[:,:H_N])
df[:,:conv] = [loaddistribution(Phash).metadata[:minimizer_converged] for Phash in df[:,:P_2]]

names(df)

for subdf in groupby(df, [:ori_mat_file, :ori_mat_rec, :neuron_type])
    println("$(ts()) File $(subdf[1,:ori_mat_file]) recording #$(subdf[1,:ori_mat_rec]) [$(subdf[1,:neuron_type]), $(extrema(subdf[:,:n_neurons])) min/max #neurons]")
    I_2 = subdf[:, :I_2]
    I_2_conv = [r[:I_2] for r in eachrow(subdf) if r[:conv]]
    N_converged = length(I_2I_2_conv)
    println("$(ts())   Found $N_converged distributions where fitting loglikelihood converged.")
    if N_converged > 0
        min_I_2, mean_I_2, max_I_2 = (minimum(I_2_conv), mean(I_2_conv), maximum(I_2_conv))
        print("$(ts())     min  I_2: "); print_with_color(0 <= min_I_2 <= 1 ? :default : :red, min_I_2); println()
        print("$(ts())     mean I_2: "); print_with_color(0 <= mean_I_2 <= 1 ? :default : :red, mean_I_2); println()
        print("$(ts())     max  I_2: "); print_with_color(0 <= max_I_2 <= 1 ? :default : :red, max_I_2); println()
    end
    # I_2 = subdf[:,:I_2]
    # print("\tmin  I_2: "); min_I_2 = minimum(skipmissing(I_2))
    # print_with_color(min_I_2 <= 0.0 ? :red : :default, min_I_2)
    # println()
    # print("\tmean I_2: "); mean_I_2 = mean(skipmissing(I_2)); mean_I_2_skip_bad = mean(I_2[0 .<= I_2 .<= 1])
    # print_with_color(mean_I_2 <= 0.0 || mean_I_2 >= 1.0 ? :red : :default, mean_I_2)
    # print("\t(cleaned) $mean_I_2_skip_bad")
    # println()
    # print("\tmax  I_2: "); max_I_2 = maximum(skipmissing(I_2))
    # print_with_color(max_I_2 >= 1.0 ? :red : :default, max_I_2)
    # println()
end
