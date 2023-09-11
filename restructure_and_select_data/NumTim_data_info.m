function stat = NumTim_data_info(stat,all_combis,select_ids,DT_runs,run)
%% Selects general info from requested ids and assigns them to stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% stat: structure containing statistical info from all the performed tests
% all_combis: struct with information on the provided DT_runs. It includes
%       information on the size, shared coordinates and their topographic
%       vectors of the compared maps. Its field sizes
%       are [amount of maps DT_run_1] (e.g. 11 for TimingAll) x [amount of
%       maps DT_run_2] (e.g., 9 for NumerosityEven) x [amount of subjects]
%       (in our case 8) x [amount of hemispheres] (in our case 2)
% select_ids: ids indicating which combinations to include in the analyses
% DT_runs: the two data runs you want to compare (options: TimingAll, TimingEven,
%       TimingOdd, NumerosityAll, NumerosityEven, NumerosityOdd)
% run: index of DT_runs
% stat_path: path where you want stat to be saved
%
% Output
% stat: structure containing statistical info from all the performed tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Find indices requested combinations
corresponding_map = strcat('Map',char(num2str(run)));
other_map = strcat('Map',char(num2str(3-run)));

% select general info corresponding to the indices
stat.data.(DT_runs{run}).subj = all_combis.subject(select_ids);
stat.data.(DT_runs{run}).map = all_combis.(corresponding_map)(select_ids);
stat.data.(DT_runs{run}).other_map = all_combis.(other_map)(select_ids);
stat.data.(DT_runs{run}).hemi= all_combis.hemi(select_ids);
