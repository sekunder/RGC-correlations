# RGC-correlations
Tools for studying information content of retinal ganglion cell activity

# Version History
I have a loose notion of "version" for this software that will eventually become better-defined.

v0.2 Changed some fields in `GrayScaleStimulus` from `Tuple` to `Vector` because Julia v0.5 doesn't have `./` defined for `Tuple` types, and I'm throwing these scripts back and forth between my laptop (running Julia 6) and a server (running Julia 5)

v0.1 First working version of STRF computation and simulation.
