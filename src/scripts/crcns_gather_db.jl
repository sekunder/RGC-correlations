# lf_name = joinpath(homedir(),"julia",basename(@__FILE__) * ".log")
help_string = """
$(basename(@__FILE__))

Systematically go through existing spike trains/rasters, probability distributions, and
STRFs that have been computed and then create "databases"; i.e. `DataFrame`s stored
in .jld files in [RGC-Correlations]/data.
"""

println(help_string)
# lf = open(lf_name, "a")
println("Initializing...")
include("../util/init.jl")
using JLD, MAT, DataFrames

# load database(s) from ../../data, or create them if they don't exist
db_dir = joinpath(@__DIR__, "../../data")
if !ispath(db_dir)
    println("Creating directory db_dir")
    mkpath(db_dir)
end

function new_strf_dataframe()
    return DataFrame(:ori_mat_file => [],
        :ori_mat_rec => [],
        :ori_mat_neuron => [],
        :ori_jld_file => [],
        :neuron_type => [],
        :hash => [],
        :frame_rate => [],
        :resolution => [],
        :pixels => [],
        :center => [],
        :covariance => [],
        :estimation_params => [])
end

function new_prob_dataframe()
    return DataFrame(:ori_mat_file => [],
        :ori_mat_rec => [],
        :ori_jld_file => [],
        :hash => [],
        :neurons => [],
        :type => [],
        :entropy => [],
        :entropy2 => [],
        :mu => [])
end

# db filenames are in CRNCS_db_X_Y where X is one of prob, strf and Y is one of real, sim
function populatestrfdataframe_and_save(input_dir, output_dir, output_filename, neuron_type)
    df = new_strf_dataframe()
    jld_files = filter(x -> endswith(x,".jld"), readdir(input_dir))
    println("$(ts()) Processing STRFs in $input_dir")
    println("$(ts()) Found these files:")
    println("$(ts()) \t$(join(jld_files,"$(ts()) \t\n"))")
    for jf in jld_files
        println("$(ts())   File $jf: ")
        STRFs = load(joinpath(CRCNS_STRF_dir,"real",jf), "STRFs")
        for (i,S) in enumerate(STRFs)
            h = hash(S)
            f = savestimulus(S; dir=joinpath(CRCNS_STRF_dir,"real"))
            mf = jf[1:11] * ".mat"
            mr = parse(jf[13:13]) #funny quirk: jf[13] is a char, hence does not parse. jf[13:13] is a string, hence can be parsed.
            push!(df, [mf mr i jf neuron_type h frame_rate(S) resolution(S) frame_size(S) missing missing missing])
            println("$(ts())     $i: $f")
        end
    end
    categorical!(df, [:ori_mat_file, :ori_jld_file, :neuron_type])
    show(head(df))
    println("$(ts()) Saving to $db_dir")
    save(joinpath(db_dir, CRCNS_db_strf_real), "db", df)
end


# First, let's get all the real STRFs
populatestrfdataframe_and_save(joinpath(CRCNS_STRF_dir,"real"), db_dir, CRCNS_db_strf_real, "mouse RGC")
# and the simulated:
populatestrfdataframe_and_save(joinpath(CRCNS_STRF_dir,"sim"), db_dir, CRNCS_db_strf_sim, "simulated RGC")
