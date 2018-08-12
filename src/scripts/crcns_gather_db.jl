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
using JLD, MAT, DataFrames, CSV

# load database(s) from ../../data, or create them if they don't exist
data_dir = normpath(joinpath(@__DIR__, "../../data"))
if !ispath(data_dir)
    println("Creating directory $data_dir")
    mkpath(data_dir)
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
function populatestrfdataframe_and_save(input_dir, db_dir, db_filename; neuron_type="mouse RGC")
    df = new_strf_dataframe()
    jld_files = filter(x -> startswith(x,"2008") && endswith(x,".jld"), readdir(input_dir))
    println("$(ts()) Processing files in $input_dir")
    println("$(ts()) Found these files:")
    println("$(ts()) \t$(join(jld_files,"\n$(ts()) \t"))")
    for jf in jld_files
        println("$(ts()) * File $jf: ")
        try
            STRFs = load(joinpath(CRCNS_STRF_dir,"real",jf), "STRFs")
            for (i,S) in enumerate(STRFs)
                try
                    h = hash(S)
                    f = savestimulus(S)
                    mf = jf[1:11] * ".mat"
                    mr = parse(jf[13:13]) #funny quirk: jf[13] is a char, hence does not parse. jf[13:13] is a string, hence can be parsed.
                    push!(df, [mf, mr, i, jf, neuron_type, h, frame_rate(S), resolution(S), frame_size(S), missing, missing, missing])
                    println("$(ts())     $i: $f")
                catch en
                    println("$(ts()) !   $i: Error: $en")
                    println("$(ts())       If STRF successfully saved, it is here: $f")
                    continue
                end
            end
        catch es
            println("$(ts()) !!! Error encountered with $jf")
            show(es)
            println()
            continue
        end
    end
    categorical!(df, [:ori_mat_file, :ori_jld_file, :neuron_type])
    println("$(ts()) $(size(df)) DataFrame created")
    show(head(df))
    println()
    println("$(ts()) Saving to $(joinpath(db_dir, db_filename*".csv"))")
    CSV.write(joinpath(db_dir,db_filename*".csv"), df)
    # save(joinpath(db_dir, db_filename), "db", df)
end


# First, let's get all the real STRFs
println("Going for real STRFs")
populatestrfdataframe_and_save(joinpath(CRCNS_STRF_dir,"real"), data_dir, CRCNS_db_strf_real; neuron_type="mouse RGC")
println()
println("Done with real STRFs!")

println("And now simulated:")
populatestrfdataframe_and_save(joinpath(CRCNS_STRF_dir,"sim"), data_dir, CRCNS_db_strf_sim; neuron_type="simulated RGC")
# and the simulated:
# populatestrfdataframe_and_save(joinpath(CRCNS_STRF_dir,"sim"), data_dir, CRNCS_db_strf_sim, "simulated RGC")
