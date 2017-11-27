
base_dir = ispath(joinpath(homedir(),"Documents","MATLAB")) ? joinpath(homedir(),"Documents","MATLAB") : homedir()
# CRCNS_dir = ispath("/Users/sekunder/Documents/MATLAB/crcns_ret-1") ? "/Users/sekunder/Documents/MATLAB/crcns_ret-1" : "/data1/homes/abk170/crcns_ret-1"
CRCNS_dir = joinpath(base_dir, "crcns_ret-1")
CRCNS_data_dir = joinpath(CRCNS_dir, "Data")
CRCNS_analysis_dir = joinpath(CRCNS_dir, "analysis")
CRCNS_STRF_dir = joinpath(CRCNS_analysis_dir, "STRF")
CRCNS_information_dir = joinpath(CRCNS_analysis_dir, "information")
CRCNS_plots_dir = joinpath(CRCNS_analysis_dir, "plots")

CRCNS_script_version = v"0.3"

Ganmor_dir = joinpath(base_dir, "ganmor-retina-data")
Ganmor_analysis_dir = joinpath(Ganmor_dir, "analysis")
Ganmor_STRF_dir = joinpath(Ganmor_analysis_dir, "STRF")
Ganmor_information_dir = joinpath(Ganmor_analysis_dir, "information")

Ganmor_script_version = v"0.1"
