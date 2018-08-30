# STATUS Aug 28: I believe I have all the strf-db-related functions working, so I should be
# able to get the others up to speed quickly. Then it's time to do the "actual work"


# functions for dealing with the "databases" I am making

@everywhere using DataFrames, CSV

# """
#     vec2str(v) = join(v, " ")
#
# Used for reading/writing CSV files (transforms vectors in `DataFrame`s)
# """
@everywhere vec2str(v) = join(v, " ") # jklol this function ends up not being used anywhere
# vec2str(::Missings.Missing) = missing

# """
#     str2vec(s) = map(parse, split(s, " "))
#
# Used for reading/writing CSV files (transforms strings in CSV files)
# """
@everywhere str2vec(s) = map(parse, split(replace(s, ['[',',',']'], "")))
@everywhere str2vec(::Missings.Missing) = missing

@everywhere function new_strf_dataframe()
    return DataFrame(
        :ori_mat_file   => [],
        :ori_mat_rec    => [],
        :ori_mat_neuron => [],
        :ori_jld_file   => [],
        :neuron_type    => [],
        :hash           => [],
        :frame_rate     => [],

        # turns out CSV.write and CSV.read have a transforms kwarg so I can make my life
        # easier here by just storing actual vectors.
        :resolution     => [],
        :pixels         => [],
        :center         => [],
        :v1             => [], # eigenvectors of covariance matrix
        :v2             => [],
        # :estimation_params => [] # I don't think this is necessary to store in db
    )
end


@everywhere function new_prob_dataframe()
    return DataFrame(
        :ori_mat_file   => [],
        :ori_mat_rec    => [],
        :ori_jld_file   => [],
        :neuron_type    => [],
        :neurons        => [],
        # due to issues with CSV, this should always be stored as an int. Use
        # int_to_index_set when performing queries and the like.

        :P_1            => [], #store the hash of the three distributions for this index set
        :P_2            => [],
        :P_N            => [],
        :H_1            => [], #store the computed entropy when available
        :H_2            => [],
        :H_N            => [],
    )
end

@everywhere function new_spikes_dataframe()
    return DataFrame(
        :ori_mat_file   => [],
        :ori_mat_rec    => [],
        :ori_jld_file   => [],
        :neuron_type    => [],
        :neurons        => [], # see comment in new_prob_dataframe; this'll be an int
        :hash           => [],
        :stimulus       => [],
    )
end

#seems to be working
@everywhere function load_strf_db(filename; dir=CRCNS_analysis_dir)
    _fn = endswith(filename, ".csv") ? filename : filename * ".csv"
    df = CSV.read(joinpath(dir, _fn),
        types=Dict(
            "ori_mat_file"      => String,
            "ori_mat_rec"       => Int,
            "ori_mat_neuron"    => Int,
            "ori_jld_file"      => Union{String,Missing},
            "neuron_type"       => String,
            "hash"              => UInt,
            "frame_rate"        => Float64,
            ),
        transforms=Dict(
            8   => str2vec,
            9   => str2vec,
            10  => str2vec,
            11  => str2vec,
            12  => str2vec,
            )
        )
    # now, a problem. if a column is all missing, the eltype of that column is set to
    # `Missing` instead of `Union{Missing, X}`. This should only be a problem for the
    # `Vector{X}` columns, which cannot be parsed as `Vector`s directly. So, lets fix that
    # manually.
    types=Dict(
        # "ori_mat_file"      => String,
        # "ori_mat_rec"       => Int,
        # "ori_mat_neuron"    => Int,
        # "ori_jld_file"      => String,
        # "neuron_type"       => String,
        # "hash"              => UInt,
        # "frame_rate"        => Float64,
        "resolution"        => Vector{Int},
        "pixels"            => Vector{Int},
        "center"            => Vector{Float64},
        "v1"                => Vector{Float64},
        "v2"                => Vector{Float64},
        # "estimation_params" => String
        )
    for (k,t) in zip(names(df), eltypes(df))
        if t == Missing
            df[k] = missings(types[string(k)], length(df[k]))
        end
    end
    return df
end

# seems to be working
@everywhere function load_prob_db(filename; dir=CRCNS_analysis_dir)
    _fn = endswith(filename, ".csv") ? filename : filename * ".csv"
    df = CSV.read(joinpath(dir, _fn),
        types=Dict(
            "ori_mat_file"  => String,
            "ori_mat_rec"   => Int,
            "ori_jld_file"  => Union{String,Missing},
            "neuron_type"   => String,
            "neurons"       => Int,
            "P_1"           => UInt,
            "P_2"           => UInt,
            "P_N"           => UInt,
            "H_1"           => Union{Float64,Missing},
            "H_2"           => Union{Float64,Missing},
            "H_N"           => Union{Float64,Missing},
            ),
        )
    # all types are specified in CSV.read, so no need to get fancy as in the strf function.
    return df
end

@everywhere function load_spikes_db(filename; dir=CRCNS_analysis_dir)
    _fn = endswith(filename, ".csv") ? filename : filename * ".csv"
    df = CSV.read(joinpath(dir, _fn),
            types=Dict(
                "ori_mat_file"  => String,
                "ori_mat_rec"   => Int,
                "ori_jld_file"  => Union{String,Missing},
                "neuron_type"   => String,
                "neurons"       => Int,
                "hash"          => UInt,
                "stimulus"      => Union{UInt, Missing}
                )
            )
    return df
end

# because of bugs with the transforms kwarg, all three of these save functions are actually
# identical. they weren't supposed to be! oh well.
@everywhere function save_strf_db(df, filename; dir=CRCNS_analysis_dir)
    _fn = endswith(filename, ".csv") ? filename : filename * ".csv"
    CSV.write(joinpath(dir, _fn), df)
end

@everywhere function save_prob_db(df, filename; dir=CRCNS_analysis_dir)
    _fn = endswith(filename, ".csv") ? filename : filename * ".csv"
    CSV.write(joinpath(dir, _fn), df)
end

@everywhere function save_spikes_db(df, filename; dir=CRCNS_analysis_dir)
    _fn = endswith(filename, ".csv") ? filename : filename * ".csv"
    CSV.write(joinpath(dir, _fn), df)
end
