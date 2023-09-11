function stat = NumTime_data_prop(stat,whichCombi,all_combis,select_ids,DT_runs,run)
%% Selects proportions from requested ids and assigns them to stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% stat: structure containing statistical info from all the performed tests
% whichCombi: 1 = compare to "All" maps; 2 = compare to most overlapping
%       map; 3 = compare to same map; 4 = compare to all overlapping maps
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
%
% Output
% stat: structure containing statistical info from all the performed tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% chance level proportion computations
% (only conceptually makes sense for whichCombi = 1 or whichCombi = 4)
all_gray_voxels = all_combis.size_total;

% To compute chance level: all responsize voxels
% for example: For numerosity maps you want to know how much overlap you get by
% chance with timing. So you compute the chance level of timing (to know
% for numerosity maps)
if whichCombi == 1
    % selection ids will get "All" map
    corresponding_to_model_1 = all_combis.size_Map1;
    corresponding_to_model_2 = all_combis.size_Map2;
elseif whichCombi == 4
    % one number per row or column, but will be in the shape of all_combis
    % when divided by all_gray_voxels
    corresponding_to_model_1 = nansum(all_combis.size_Map1(1:end-1,:,:,:),1);
    corresponding_to_model_2 = nansum(all_combis.size_Map2(:,1:end-1,:,:),2);
end

% works for whichCombi 1 & 4 as the ids indicate the "All" map (which is
% also where general info was selected) and chance proportion is everywhere
% the same for each hemisphere
if whichCombi == 1 || whichCombi == 4
    eval(['chance_proportions = corresponding_to_model_', char(num2str(3-run)),'./all_gray_voxels;'])
    stat.data.(DT_runs{run}).chance_proportion = chance_proportions(select_ids);
end

%% size maps and overlap
% size of the selected maps

eval(['stat.data.(DT_runs{run}).size = all_combis.size_Map',char(num2str(run)),'(select_ids);'])

% proportions of overlap
amount_overlap = arrayfun(@(x) length(x{:}), all_combis.shared_coord_ids);

if whichCombi == 4
    % could be done way nicer, but now I know it works correctly
    for datapoint = 1:length(stat.data.(DT_runs{run}).subj)
        eval(['total_overlap(datapoint) = sum(amount_overlap(all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & all_combis.Map',char(num2str(run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(all_combis.Map',char(num2str(3-run)),', "All")));'])
    end
else
    total_overlap = amount_overlap(select_ids)';
end

stat.data.(DT_runs{run}).prop = total_overlap'./stat.data.(DT_runs{run}).size;