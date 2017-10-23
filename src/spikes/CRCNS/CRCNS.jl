# Utility functions for dealing with CRCNS data specifically.

using MAT

"""
    get_spikes_from_file(filename, rec_idx)

Convenience function for extracting spike times from the given matlab file from
the CRCNS data.

"""
function CRCNS_get_spikes_from_file(filename, rec_idx::Int)
    vars = matread(filename)

    spike_trains = vars["spikes"][:,rec_idx]
    return SpikeTrains(spike_trains;
        stim_start=vars["stimulus"]["onset"][rec_idx],
        stim_finish=vars["stimulus"]["onset"][rec_idx] + vars["stimulus"]["frame"][rec_idx] * vars["stimulus"]["Nframes"][rec_idx],
        filename=filename)
end
