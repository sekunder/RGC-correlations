
@everywhere base_dir = ispath(joinpath(homedir(),"Documents","MATLAB")) ? joinpath(homedir(),"Documents","MATLAB") : homedir()
# CRCNS_dir = ispath("/Users/sekunder/Documents/MATLAB/crcns_ret-1") ? "/Users/sekunder/Documents/MATLAB/crcns_ret-1" : "/data1/homes/abk170/crcns_ret-1"
@everywhere CRCNS_dir = joinpath(base_dir, "crcns_ret-1")
@everywhere CRCNS_data_dir = joinpath(CRCNS_dir, "Data")
@everywhere CRCNS_analysis_dir = joinpath(CRCNS_dir, "analysis")
@everywhere CRCNS_STRF_dir = joinpath(CRCNS_analysis_dir, "STRF")
@everywhere CRCNS_information_dir = joinpath(CRCNS_analysis_dir, "information")
@everywhere CRCNS_plots_dir = joinpath(CRCNS_analysis_dir, "plots")

@everywhere CRCNS_script_version = v"0.4"

@everywhere Ganmor_dir = joinpath(base_dir, "ganmor-retina-data")
@everywhere Ganmor_analysis_dir = joinpath(Ganmor_dir, "analysis")
@everywhere Ganmor_STRF_dir = joinpath(Ganmor_analysis_dir, "STRF")
@everywhere Ganmor_information_dir = joinpath(Ganmor_analysis_dir, "information")

@everywhere Ganmor_script_version = v"0.1"
