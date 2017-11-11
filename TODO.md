# Project Tasks
A big-picture view of features to implement and scripts to write; concrete TODO items if appropriate.

# STRF Clean up

Some of the neurons from the CRCNS recordings seem to be noise. That is, the
computed STRF does not have a a clear "center", meaning they don't play into the
geometric argument we are putting forth. So, I need a way to filter those out.

# STRF geometric analysis

In that same vein, I need to compute the geometric properties of the STRFs I'm
computing -- RF center, things like that.
[ ] function: (grayscale stimulus) -> geometric info (see STRF_image_stats in [convexneuron]src/spikes/CRCNS_STRF.jl)

# Information-theoretic analysis

I need scripts to run the information-theoretic analysis of the real and
simulated data. Most of the tools for this are already written, I just need to
glue them all together.
[ ] Convenience func 1: (spike trains, index set, output file) -> write P_1, P_2, P_N to output file, computing entropies when possible
[ ] Convenience func 2: (CRCNS file name, sample size, num trials) -> loop func 1 several times, using both real and sim data
[ ] Script: loop through all CRCNS files; check if simulated data exists; run func 2 over a range of sample sizes

# Plots!

I am following a general trend of scripts that perform one task (take data,
perform one specific manipulation, dump it to file) so this will involve writing
scripts that read these giant collections of data and produce plots.

# Ganmoor Data

I need to write a function similar to `CRCNS_Stimulus` that will produce a
`GrayScaleStimulus` object from the Ganmoor data. This shouldn't be too bad.
Likewise, a function to load the spikes... but perhaps before doing this I
should do some of the miscellaneous tasks in the General Coding Stuff section
(Ganmoor data includes "spike quality" information, which I would like to
include in "hidden metadata")

# General Coding Stuff

[ ] Think about what I need from a "version" system.
[x] Implement a "hidden metadata" system. The problem I want to solve is: I might want to include information in the `metadata` field which could be, say, a vector of length > 1000. I don't want to print that value in `show_metadata`, but maybe just the key.
    Possible solution: Make a `:hidden` key for metadata; `hide_metadata!(A, key)` just adds `key` to `A.metadata[:hidden]`, then `show_metadata` checks if there's any keys to skip.
[ ] Move everything from `/src` to top level, maybe?
[x] `probability.jl`: find func defs of the form `func(args) = _func(args)` and change them
[x] `GrayScaleStimulus.jl` add information about resolution to `show`
[ ] `stimulus.jl` move default_CRCNS_dir definition to new file `/util/init.jl`
