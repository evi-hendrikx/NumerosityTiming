function stat = NumTim_prepare_whichCombi(all_combis, whichCombi, topo_measurement_per_map,minAmountCorrelation,type_angle,select_ids,DT_runs,stat_path,selection_all_combis,selection_DT_runs)
%% Selects datapoints for requested combination and assigns them to stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% all_combis: struct with information on the provided DT_runs. It includes
%       information on the size, shared coordinates and their topographic
%       vectors of the compared maps. Its field sizes
%       are [amount of maps DT_run_1] (e.g. 11 for TimingAll) x [amount of
%       maps DT_run_2] (e.g., 9 for NumerosityEven) x [amount of subjects]
%       (in our case 8) x [amount of hemispheres] (in our case 2)
% whichCombi: 1 = compare to "All" maps; 2 = compare to most overlapping
%       map; 3 = compare to same map; 4 = compare to all overlapping maps;
%       5 = use voxel selection from other data structure
% topo_measurement_per_map: calculate mean or median angle between vectors
% minAmountCorrelation: minimum amount of voxels a map should have in order
%       to calculate the correlation coefficients (and topo angles)
% type_angle: angles between topographic vectors can be computed taking
%       into account that the hemispheres are mirror images of each other
%       (medial-lateral), or using the exact same method for both
%       hemispheres (left-right)
% select_ids: ids indicating which combinations to include in the analyses
% DT_runs: the two data runs you want to compare (options: TimingAll, TimingEven,
%       TimingOdd, NumerosityAll, NumerosityEven, NumerosityOdd)
% stat_path: path where you want stat to be saved
%
% Output
% stat: structure containing statistical info from all the performed tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% select appropriate data for whichCombi
stat = {};

% you want the selection of the same variable (e.g., comparing NumEven and
% NumOdd, you want to use all overlapping maps for the numerosity maps, not
% the timing maps)
if whichCombi == 5
    if contains(DT_runs{1},'Numerosity') == 1
        selection_run = find(contains(selection_DT_runs,'Numerosity'));
    else
        selection_run = find(contains(selection_DT_runs,'Timing'));
    end
end

for run = 1:length(DT_runs)
    stat = NumTim_data_info(stat,all_combis,select_ids{run},DT_runs,run);

    if whichCombi == 5
        stat = NumTime_data_corr(stat,whichCombi,all_combis,select_ids{run},DT_runs,run,minAmountCorrelation,selection_all_combis,selection_DT_runs,selection_run);
        stat = NumTime_data_topo(stat,whichCombi,all_combis,select_ids{run},DT_runs,run,minAmountCorrelation,topo_measurement_per_map,type_angle,stat_path,selection_all_combis,selection_DT_runs,selection_run);

    else
        stat = NumTime_data_prop(stat,whichCombi,all_combis,select_ids{run},DT_runs,run);
        stat = NumTime_data_corr(stat,whichCombi,all_combis,select_ids{run},DT_runs,run,minAmountCorrelation);
        stat = NumTime_data_topo(stat,whichCombi,all_combis,select_ids{run},DT_runs,run,minAmountCorrelation,topo_measurement_per_map,type_angle,stat_path);

    end
end

if whichCombi ~= 5 && whichCombi ~= 1
    stat_path = char(stat_path);
    tbl_path = insertAfter(stat_path, 'stat','_topoTbl');
    tbl_path = [tbl_path(1:end-3), 'csv'];
    if (contains(DT_runs{1},'Numerosity') == 1 && contains(DT_runs{2},'Numerosity') == 1) || (contains(DT_runs{1},'Timing') == 1 && contains(DT_runs{2},'Timing') == 1)
        create_table_for_R(stat, tbl_path,1);
    else
        create_table_for_R(stat, tbl_path,0);
    end
end


