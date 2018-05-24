# Utility functions for dealing with CRCNS data specifically.

using MAT

"""
    CRCNS_get_spikes_from_file(filename, rec_idx)

Convenience function for extracting spike times from the given matlab file from
the CRCNS data.

"""
@everywhere function CRCNS_get_spikes_from_file(filename, rec_idx::Int)
    vars = matread(filename)

    spike_trains = vars["spikes"][:,rec_idx]
    spikes = [st - vars["stimulus"]["onset"][rec_idx] for st in spike_trains]
    return SpikeTrains(spikes;
        stim_start=0.0,
        stim_finish=vars["stimulus"]["frame"][rec_idx] * vars["stimulus"]["Nframes"][rec_idx],
        filename=filename)
end
