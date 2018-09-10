using PyPlot

include("../util/init.jl")

# lfname = joinpath(homedir(), "julia", basename(@__FILE__) * ".log")
# lf = open(lfname, "a")

raster_bin_size = 0.020

"""
    selectcols_rename(df; keep=names(df), drop=[], name=Dict())

Returns a copy of `df`, returning only the columns in
`setdiff(union(keep, keys(name)), drop)`, and renaming columns according to the pairs in
`name`

"""
function selectcols_rename(df; keep=names(df), drop=Symbol[], name=Dict{:Symbol,:Symbol}())
    cs = setdiff(union(keep, keys(name)), drop)
    _df = copy(df[:, cs])
    rename!(_df, name)
end

################################################################################
### Spikes database stuff
################################################################################

println("$(ts()) Loading spikes databases")
dfspikes_real = load_spikes_db(CRCNS_db_spikes_real)
dfspikes_sim = load_spikes_db(CRCNS_db_spikes_sim)

spikes_real = selectcols_rename(dfspikes_real,
    drop=[:ori_jld_file, :neurons, :stimulus, :n_neurons],
    name=Dict(:hash=>:hash_real))

spikes_sim = selectcols_rename(dfspikes_sim,
    drop=[:ori_jld_file, :neuron_type, :neurons, :stimulus],
    name=Dict(:hash=>:hash_sim))

spikes_db = join(spikes_real, spikes_sim, on=[:ori_mat_file, :ori_mat_rec])

################################################################################
### STRF database stuff
################################################################################
println("$(ts()) Loading STRF databases")
dfstrf_real = load_strf_db(CRCNS_db_strf_real)
dfstrf_sim = load_strf_db(CRCNS_db_strf_sim)

# names(dfstrf_real)

strfs_real = selectcols_rename(dfstrf_real,
    drop=[:ori_jld_file, :frame_rate, :resolution, :pixels, :mm_per_px],
    name=Dict(:hash=>:hash_real, :center=>:center_real, :v1=>:v1_real, :v2=>:v2_real))

strfs_sim = selectcols_rename(dfstrf_sim,
    drop=[:ori_jld_file, :neuron_type],
    name=Dict(:hash=>:hash_sim, :center=>:center_sim, :v1=>:v1_sim, :v2=>:v2_sim))

strf_db = join(strfs_real, strfs_sim, on=[:ori_mat_file,:ori_mat_rec,:ori_mat_neuron])

strf_db[:c_diff] = missings(Float64, size(strf_db,1)); strf_db[:s_diff] = missings(Float64, size(strf_db,1));

println("$(ts()) Computing STRF differences")
for r in eachrow(strf_db)
    rf_real = loadstimulus(r[:hash_real])
    rf_sim = loadstimulus(r[:hash_sim])
    r[:s_diff] = norm(matrix_form(rf_real) - matrix_form(rf_sim))
    r[:c_diff] = (ismissing(r[:center_real]) || ismissing(r[:center_sim])) ? missing : norm(r[:center_real] - r[:center_sim])
end

# sdiffs = zeros(0); cdiffs=zeros(0);
# for r in eachrow(strf_db[:, [:s_diff, :c_diff]])
#     if !(ismissing(r[:s_diff]) || ismissing(r[:c_diff]))
#         push!(sdiffs, r[:s_diff])
#         push!(cdiffs, r[:c_diff])
#     end
# end
# scatter(sdiffs, cdiffs)

# rf_diffs = rename!(by(strf_db,
#     [:ori_mat_file, :ori_mat_rec],
#     df -> [mean(skipmissing(df[:s_diff])) mean(skipmissing(df[:c_diff]))]
#     ), :x1=>:mean_sdiff, :x2=>:mean_cdiff)


prob_db = new_prob_dataframe()

for subdf in groupby(strf_db, [:ori_mat_file, :ori_mat_rec])
    mf = subdf[1,:ori_mat_file]; mr = subdf[1,:ori_mat_rec]; N_neurons = size(subdf,1)
    println("$(ts()) Looking at data file $mf, recording $mr. Found $N_neurons neurons.")
    # println(describe(subdf[:,[:c_diff,:s_diff]]))
    # println("var_cdiff: $(var(skipmissing(subdf[:c_diff])))")
    # println("var_sdiff: $(var(skipmissing(subdf[:s_diff])))")
    # println()

    if size(subdf, 1) >= 10
        _df = DataFrame(subdf[:,:])
        # if there's at least 10 neurons to work with, let's pull up the spikes for that
        # file/rec, and start fitting things to the 10 "best" neurons, then 11, and so on up
        # to 20.
        sp_idx = find((spikes_db[:ori_mat_file] .== mf) .& (spikes_db[:ori_mat_rec] .== mr))[1]
        # println("using this row from spikes_db:")
        # println(spikes_db[sp_idx, :])
        println("$(ts()) Preparing to fit just a whole mess of probability distributions")
        X_real = raster(loadspikes(spikes_db[sp_idx, :hash_real]), raster_bin_size)
        X_sim = raster(loadspikes(spikes_db[sp_idx, :hash_sim]), raster_bin_size)
        # println("Since ")

        sort!(_df, [:s_diff])
        for N = 10:min(20, N_neurons)
            neurons = sort(subdf[1:N, :ori_mat_neuron])
            println("$(ts())   Using the $N neurons with lowest STRF difference")
            println("$(ts())   neurons: $neurons")
            # println("Testing: Would try to fit on neurons $neurons")
            for (nt,X) in zip([_df[1,:neuron_type],"simulated"], [X_real, X_sim])
                P_1 = first_order_model(X, neurons); savedistribution(P_1)
                println("$(ts())     $nt\tP_1: $(hash(P_1))")
                P_2 = second_order_model(X, neurons); savedistribution(P_2)
                println("$(ts())     $nt\tP_2: $(hash(P_2))")
                P_N = data_model(X, neurons); savedistribution(P_N)
                println("$(ts())     $nt\tP_N: $(hash(P_N))")
                H_1 = entropy2(P_1)
                H_2 = entropy2(P_2)
                H_N = entropy2(P_N)
                push!(prob_db,
                    [mf,
                    mr,
                    nt,
                    N,
                    index_set_to_int(neurons),
                    hash(P_1),
                    hash(P_2),
                    hash(P_N),
                    H_1,
                    H_2,
                    H_N])
            end
        end
    end
end

println("$(ts()) Saving probability distributions database")
save_prob_db(prob_db, CRCNS_db_prob)

################################################################################
### Experimenting with stuff
################################################################################

# for r in eachrow(spikes_db)
#     mf = r[:ori_mat_file]; mr = r[:ori_mat_rec];
#     # i could use Query.jl to do some querying here but whatever.
#     # I want all the rows of strf_db where the file is mf and the recording is mr
#     rows = find((strf_db[:ori_mat_file] .== mf) & (strf_db[:ori_mat_rec] .== mr))
#     strfs = strf_db[rows, :]
#     sort!(strfs, [:s_diff])
# end
