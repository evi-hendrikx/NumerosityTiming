%% Runs the entire pipeline for NumTime
%% necessary variables
clear all
close all

general_folder = '[gen_path]';
addpath(general_folder)
addpath(strcat(general_folder, '/scripts'));
addpath(strcat(general_folder, '/scripts/load_data'));
addpath(strcat(general_folder, '/scripts/restructure_and_select_data'));
addpath(strcat(general_folder, '/scripts/statistics'));
addpath(strcat(general_folder, '/scripts/PyColormap4Matlab'));
addpath(strcat(general_folder, '/scripts/circstat-matlab-master'));
addpath(strcat(general_folder, '/scripts/repeated_measures_comparisons'));

paths{1}=[anonymized]
paths{2}=[anonymized]
paths{3}=[anonymized]
paths{4}=[anonymized]
paths{5}=[anonymized]
paths{6}=[anonymized]
paths{7}=[anonymized]
paths{8}=[anonymized]

new_subjNames = {'S6','S5','S4','S7','S8','S3','S1','S2'};
data_path = strcat(general_folder,'/data');
save_path = strcat(general_folder,'/results');

% for selecting data
minVEnum=0.3;
minVEtime=0.2;

%% Load in All, Even, and Odd scans for all ROIs and all models into a struct
NumTim_path = strcat(data_path, '/NumTim_data_minVEtime=',string(minVEtime),'_minVEnum=',string(minVEnum),'.mat');
NumTim_data = NumTim_load_data(NumTim_path, new_subjNames, paths, minVEtime, minVEnum);

if exist(strcat(data_path, '/time_series.mat'), 'file') ~= 2 % Only needs to be done if you don't already have time_series_
    NumTim_time_series(data_path, new_subjNames, paths)
end

%% Compute topographic info and restructure data
% for determining compared pairs// Keep alphabetical order
DT_run_1 = 'NumerosityAll';
DT_run_2 = 'TimingAll';
DT_runs = {DT_run_1, DT_run_2};

all_combi_path = strcat(save_path, '/all_combi_',DT_run_1,'_',DT_run_2,'_minVEtime=',char(string(minVEtime)),'_minVEnum=',char(string(minVEnum)),'.mat');
topo_path = strcat(save_path, '/topo_',DT_run_1,'_',DT_run_2, '_minVEtime=',char(string(minVEtime)),'_minVEnum=',char(string(minVEnum)),'.mat');

% compute relevant aspects for all possible combinations
all_combis = NumTim_restructure_info(all_combi_path, topo_path, NumTim_data, paths, DT_runs);
%% Compute statistics for relevant combinations

% specify types of comparisons
% 1 = all overlap; 
% 2 = most overlapping map; 
% 3 = same map (e.g., NTO even vs. NTO Odd); 
% 4 = all overlapping maps; 
% 5 = use voxel selection from other data structure (both use whichCombi = 4) 
% (e.g., NTO Even vs. Numerosity maps odd, but for overlapping voxels of NTO Even vs. Timing maps odd)
% --> latter not possible for all scripts (used to compare correlations)
whichCombi = 4; 

% get selection ids to use at all_combis structure
select_ids = NumTim_get_select_ids(all_combis, whichCombi,DT_runs);

% get the information about the other selection you're interested in
if whichCombi == 5
    selection_DT_run_1 = 'NumerosityEven';
    selection_DT_run_2 = 'TimingOdd';

    selection_DT_runs = {selection_DT_run_1, selection_DT_run_2};
    selection_all_combi_path = strcat(save_path, '/all_combi_',selection_DT_run_1,'_',selection_DT_run_2,'_minVEtime=',char(string(minVEtime)),'_minVEnum=',char(string(minVEnum)),'.mat');
    selection_topo_path = strcat(save_path, '/topo_',selection_DT_run_1,'_',selection_DT_run_2, '_minVEtime=',char(string(minVEtime)),'_minVEnum=',char(string(minVEnum)),'.mat');

    selection_all_combis = NumTim_restructure_info(selection_all_combi_path, selection_topo_path, NumTim_data, paths, selection_DT_runs);
end

save_path = [general_folder,'/results/whichCombi=', char(num2str(whichCombi))];

%% select correct ids and compute stats
minAmountCorrelation = 10;
type_angle_options = {'medial-lateral','left-right'};
type_angle = type_angle_options{1};
topo_measurement_per_map = "mean"; % "median" or "mean"

if whichCombi == 5
    stat_path = strcat(save_path, '/stat_',DT_run_1,'_',DT_run_2,'_minVEtime=',string(minVEtime),'_minVEnum=',string(minVEnum),'_whichCombi=',char(string(whichCombi)),'_selection=',selection_DT_run_1,'_',selection_DT_run_2,'_topo_measure=',char(topo_measurement_per_map),'.mat');
    stat = NumTim_prepare_whichCombi(all_combis, whichCombi,topo_measurement_per_map,minAmountCorrelation,type_angle,select_ids,DT_runs,stat_path,selection_all_combis,selection_DT_runs);
    save(stat_path, 'stat')
else
    stat_path = strcat(save_path, '/stat_',DT_run_1,'_',DT_run_2,'_minVEtime=',string(minVEtime),'_minVEnum=',string(minVEnum),'_whichCombi=',char(string(whichCombi)),'_topo_measure=',char(topo_measurement_per_map),'.mat');
    stat = NumTim_prepare_whichCombi(all_combis, whichCombi,topo_measurement_per_map,minAmountCorrelation,type_angle,select_ids,DT_runs,stat_path);
    stat = NumTim_chance_anova_stats(stat_path, whichCombi,stat,DT_runs);

    if whichCombi ==4
        topo_size
    end

    %% Make vizualizations
    fig_path = strcat(save_path, '/','minVEtime=',string(minVEtime),'_minVEnum=',string(minVEnum),'_whichCombi=',char(string(whichCombi)));
    NumTim_plots_prop_corr(fig_path, NumTim_data, stat, whichCombi, DT_runs)
    NumTim_plots_topo(fig_path, NumTim_data, stat, DT_runs)

end


%% Compare strength of correlations between num-time and self
fig_path = strcat(save_path, '/','compare_rm_proportion_minVEtime=',string(minVEtime),'_minVEnum=',string(minVEnum),'_whichCombi=',char(string(whichCombi)));
NumTim_compare_proportions(fig_path)

fig_path = strcat(save_path, '/','compare_rm_correlation_minVEtime=',string(minVEtime),'_minVEnum=',string(minVEnum),'_whichCombi=',char(string(whichCombi)));
NumTim_compare_correlations(fig_path,whichCombi)

fig_path = strcat(save_path, '/fig_rmCompare_topo_meanRuns=1_Halves_minVEtime=',string(minVEtime),'_minVEnum=',string(minVEnum),'_whichCombi=',num2str(whichCombi),'_topo_measure=',topo_measurement_per_map);
NumTim_compare_topo(fig_path, NumTim_data, whichCombi)

