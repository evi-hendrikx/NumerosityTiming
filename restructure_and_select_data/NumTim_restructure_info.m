function all_combis = NumTim_restructure_info(all_combi_path, topo_path, NumTim_data, paths, DT_runs)
%% Makes all possible comparisons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% all_combi_path: where you want your all_combis structure to be saved
% topo_path: where you want the structure with topographic vectors to be
% saved
% NumTim_data: made with the NumTim_load_data script. Contains all values
%       that you could need for all following analyses
% paths: where the information about each participant is stored (should be
%       as long as and corresponding with new_subjNames)
% DT_runs: the two data runs you want to compare (options: TimingAll, TimingEven,
%       TimingOdd, NumerosityAll, NumerosityEven, NumerosityOdd)
%
%
% Output
% all_combis: struct with information on the provided DT_runs. It includes
%       information on the size, shared coordinates and their topographic
%       vectors of the compared maps. Its field sizes
%       are [amount of maps DT_run_1] (e.g. 11 for TimingAll) x [amount of
%       maps DT_run_2] (e.g., 9 for NumerosityEven) x [amount of subjects]
%       (in our case 8) x [amount of hemispheres] (in our case 2)
%
% also creates and saves
% topo: struct with normal vectors and parameter vectors for each map in
%       the requested runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check if it already exists
if (exist(all_combi_path, 'file') == 2)
    load(all_combi_path)
    return
end

%% general info and figure out what to compare
[subjNames,Hemispheres,Run_1, Run_2] = NumTim_get_info(NumTim_data,DT_runs);

%% prepare "all combis" structure to fill
topo = {};

% general info
all_combis = {};
all_combis.subject = strings(length(Run_1.Names),length(Run_2.Names),length(subjNames), length(Hemispheres));
all_combis.hemi = all_combis.subject; all_combis.Map1 = all_combis.subject; all_combis.Map2 = all_combis.subject;

% sizes
all_combis.shared_coord_ids = cell(length(Run_1.Names),length(Run_2.Names),length(subjNames), length(Hemispheres));
all_combis.size_Map1 = NaN(length(Run_1.Names),length(Run_2.Names),length(subjNames), length(Hemispheres));
all_combis.size_Map2 = all_combis.size_Map1;
all_combis.coords_Map1 = all_combis.shared_coord_ids;
all_combis.coords_Map2 = all_combis.shared_coord_ids;

% values
all_combis.preferred_values.Run_1.x0 = cell(length(Run_1.Names),length(Run_2.Names),length(subjNames), length(Hemispheres));
all_combis.preferred_values.Run_2.x0 = all_combis.preferred_values.Run_1.x0;

% additional if include period (for timing)
if Run_1.Type == "Timing"
    all_combis.preferred_values.Run_1.y0 = all_combis.preferred_values.Run_1.x0;
end
if Run_2.Type == "Timing"
    all_combis.preferred_values.Run_2.y0 = all_combis.preferred_values.Run_1.x0;
end

%topo: normals, vectors & angles
all_combis.shared_topo_vectors.normals = all_combis.preferred_values.Run_1.x0;
all_combis.shared_topo_vectors.cross.Run_1.x0 = all_combis.preferred_values.Run_1.x0;
all_combis.shared_topo_vectors.cross.Run_2.x0 = all_combis.preferred_values.Run_1.x0;

% additional if include period (for timing)
if Run_1.Type == "Timing"
    all_combis.shared_topo_vectors.cross.Run_1.y0 = all_combis.preferred_values.Run_1.x0;
end
if Run_2.Type == "Timing"
    all_combis.shared_topo_vectors.cross.Run_2.y0 = all_combis.preferred_values.Run_1.x0;
end

%% compute shared coordinates and topographic vectors
for subj= 1:length(subjNames)
    for Hemisphere=1:length(Hemispheres)

        % so topology is only done once for each map
        do_run_2 = 1;

        disp("Start topo"); disp(subjNames{subj}); disp(Hemispheres{Hemisphere})

        for Map1 = 1:length(Run_1.Names)

            % check whether map exists
            if ~isfield(NumTim_data.(subjNames{subj}),Run_1.Names{Map1}) || ~isfield(NumTim_data.(subjNames{subj}).(Run_1.Names{Map1}),Hemispheres{Hemisphere})
                continue
            end

            % so topology is only done once for each map
            do_run_1 = 1;
            above_thresh_run_1_vol = NumTim_data.(subjNames{subj}).(Run_1.Names{Map1}).(Hemispheres{Hemisphere}).(Run_1.Type).(Run_1.DTrun).VOLUMEIndicesAboveThreshold;
            above_thresh_run_1 = NumTim_data.(subjNames{subj}).(Run_1.Names{Map1}).(Hemispheres{Hemisphere}).(Run_1.Type).(Run_1.DTrun).IndicesAboveThreshold;

            for Map2 = 1:length(Run_2.Names)

                % check whether map exists
                if ~isfield(NumTim_data.(subjNames{subj}), Run_2.Names{Map2}) || ~isfield(NumTim_data.(subjNames{subj}).(Run_2.Names{Map2}),Hemispheres{Hemisphere})
                    continue
                end

                %% fill structure with general information from loop
                all_combis.Map1(Map1,Map2,subj,Hemisphere) = string(Run_1.Names{Map1});
                all_combis.Map2(Map1,Map2,subj,Hemisphere) = string(Run_2.Names{Map2});
                all_combis.subject(Map1,Map2,subj,Hemisphere) = string(subjNames{subj});
                all_combis.hemi(Map1,Map2,subj,Hemisphere) = string(Hemispheres{Hemisphere});

                %% which voxels are above threshold & which are shared
                above_thresh_run_2_vol = NumTim_data.(subjNames{subj}).(Run_2.Names{Map2}).(Hemispheres{Hemisphere}).(Run_2.Type).(Run_2.DTrun).VOLUMEIndicesAboveThreshold;
                above_thresh_run_2 = NumTim_data.(subjNames{subj}).(Run_2.Names{Map2}).(Hemispheres{Hemisphere}).(Run_2.Type).(Run_2.DTrun).IndicesAboveThreshold;

                [shared_coord_ids, IntersectRun1, IntersectRun2]=intersectCols(above_thresh_run_1_vol, above_thresh_run_2_vol);
                all_combis.shared_coord_ids(Map1,Map2,subj,Hemisphere) = {shared_coord_ids};

                if ~isempty(IntersectRun1)
                    all_combis.preferred_values.Run_1.x0(Map1,Map2,subj,Hemisphere) = {NumTim_data.(subjNames{subj}).(Run_1.Names{Map1}).(Hemispheres{Hemisphere}).(Run_1.Type).(Run_1.DTrun).x0s(above_thresh_run_1(IntersectRun1))};
                    all_combis.preferred_values.Run_2.x0(Map1,Map2,subj,Hemisphere) = {NumTim_data.(subjNames{subj}).(Run_2.Names{Map2}).(Hemispheres{Hemisphere}).(Run_2.Type).(Run_2.DTrun).x0s(above_thresh_run_2(IntersectRun2))};
                    if Run_1.Type == "Timing"
                        all_combis.preferred_values.Run_1.y0(Map1,Map2,subj,Hemisphere) = {NumTim_data.(subjNames{subj}).(Run_1.Names{Map1}).(Hemispheres{Hemisphere}).(Run_1.Type).(Run_1.DTrun).y0s(above_thresh_run_1(IntersectRun1))};
                    end
                    if Run_2.Type == "Timing"
                        all_combis.preferred_values.Run_2.y0(Map1,Map2,subj,Hemisphere) = {NumTim_data.(subjNames{subj}).(Run_2.Names{Map2}).(Hemispheres{Hemisphere}).(Run_2.Type).(Run_2.DTrun).y0s(above_thresh_run_2(IntersectRun2))};
                    end
                end

                all_combis.size_Map1(Map1,Map2,subj,Hemisphere)=length(above_thresh_run_1_vol);
                all_combis.size_Map2(Map1,Map2,subj,Hemisphere)=length(above_thresh_run_2_vol);
                all_combis.coords_Map1(Map1,Map2,subj,Hemisphere)={above_thresh_run_1_vol};
                all_combis.coords_Map2(Map1,Map2,subj,Hemisphere)={above_thresh_run_2_vol};
                all_combis.size_total(Map1,Map2,subj,Hemisphere)=length(NumTim_data.(subjNames{subj}).GrayAll.(Hemispheres{Hemisphere}).(Run_1.Type).(Run_1.DTrun).roiIndices);

                %% get topology for all of the maps separately
                % Not doing this for the "All" maps because it takes FOREVERRRR

                if (do_run_1 == 1 && Map1 < length(Run_1.Names)) || do_run_2 == 1 && Map2 < length(Run_2.Names)

                    % doing this in a for-loop because the first time both do_runs are 1
                    run_nrs = [];
                    if do_run_1 == 1 && Map1 < length(Run_1.Names)
                        run_nrs = [run_nrs,'1'];
                    end
                    if do_run_2 == 1 && Map2 < length(Run_2.Names)
                        run_nrs = [run_nrs,'2'];
                    end

                    for run_nr_id = 1:length(run_nrs)
                        run_nr = run_nrs(run_nr_id);

                        eval(strcat('disp(Run_',run_nr,'.Names{Map',run_nr,'})'))
                        eval(strcat('topo_vector = NumTim_get_topo(paths{subj}, NumTim_data, subjNames{subj}, Hemispheres{Hemisphere}, Run_',run_nr,'.Names{Map', run_nr,'}, Run_',run_nr,'.DTrun);'))

                        % assign result to topo structure
                        eval(strcat('topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_',run_nr,'.Names{Map',run_nr,'}).included_coord_ids = topo_vector.included_coord_ids;'))
                        eval(strcat('topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_',run_nr,'.Names{Map',run_nr,'}).norm_surf = topo_vector.norm_surf;'))

                        if eval(strcat('Run_',run_nr,'.Type == "Timing"'))
                            mapTypes = {'durationMap','periodMap'};
                        else
                            mapTypes = {'numerosityMap'};
                            mT = 1;
                            eval(strcat('topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_',run_nr,'.Names{Map',run_nr,'}).(Run_',run_nr,'.DTrun).(mapTypes{mT}).mean_vector_exp  = topo_vector.(mapTypes{mT}).mean_vector_exp;'))
                            eval(strcat('topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_',run_nr,'.Names{Map',run_nr,'}).(Run_',run_nr,'.DTrun).(mapTypes{mT}).cross_mean_vector_exp = topo_vector.(mapTypes{mT}).cross_mean_vector_exp;'))
                            eval(strcat('topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_',run_nr,'.Names{Map',run_nr,'}).(Run_',run_nr,'.DTrun).(mapTypes{mT}).cross_cross_mean_vector_exp = topo_vector.(mapTypes{mT}).cross_cross_mean_vector_exp;'))
                        end

                        for mT = 1:length(mapTypes)
                            eval(strcat('topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_',run_nr,'.Names{Map',run_nr,'}).(Run_',run_nr,'.DTrun).(mapTypes{mT}).mean_vector  = topo_vector.(mapTypes{mT}).mean_vector;'))
                            eval(strcat('topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_',run_nr,'.Names{Map',run_nr,'}).(Run_',run_nr,'.DTrun).(mapTypes{mT}).cross_mean_vector = topo_vector.(mapTypes{mT}).cross_mean_vector;'))
                            eval(strcat('topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_',run_nr,'.Names{Map',run_nr,'}).(Run_',run_nr,'.DTrun).(mapTypes{mT}).cross_cross_mean_vector = topo_vector.(mapTypes{mT}).cross_cross_mean_vector;'))
                        end
                    end
                end
                save(topo_path, 'topo')

                %% get topographic vectors of overlapping coordinates
                % (not for "All" maps)
                if ~isempty(all_combis.shared_coord_ids(Map1,Map2,subj,Hemisphere)) && Map1 < length(Run_1.Names) && Map2 < length(Run_2.Names)

                    % shared coords may be in another order than topo coordinates
                    all_coord_ids_map_1 = topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_1.Names{Map1}).included_coord_ids;
                    all_coord_ids_map_2 = topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_2.Names{Map2}).included_coord_ids;
                    all_shared_coord_ids = all_combis.shared_coord_ids{Map1,Map2,subj,Hemisphere};
                    [~,~,get_vec_ids_map_1] = intersectCols(all_shared_coord_ids,all_coord_ids_map_1);
                    [~,~,get_vec_ids_map_2] = intersectCols(all_shared_coord_ids,all_coord_ids_map_2);
                    if length(all_shared_coord_ids) ~= length(get_vec_ids_map_1) || length(all_shared_coord_ids) ~= length(get_vec_ids_map_2)
                        error('topo data and proportion data have different shared coordinates')
                    end

                    % vectors for timing and numerosity have different information
                    for run = 1:2
                        if eval(strcat('Run_',char(num2str(run)),'.Type == "Timing"'))
                            eval(strcat('x_param_',char(num2str(run)),char("= 'durationMap';")))
                            eval(strcat('y_param_',char(num2str(run)),char("= 'periodMap';")))
                            eval(strcat('variable_',char(num2str(run)),char("= 'cross_mean_vector';")))
                        else
                            eval(strcat('x_param_',char(num2str(run)),char("= 'numerosityMap';")))
                            eval(strcat('variable_',char(num2str(run)),char("= 'cross_mean_vector_exp';")))
                        end
                    end

                    % select info and assign to all_combis
                    all_combis.shared_topo_vectors.cross.Run_1.x0{Map1,Map2,subj,Hemisphere} = topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_1.Names{Map1}).(Run_1.DTrun).(x_param_1).(variable_1)(:,get_vec_ids_map_1);
                    all_combis.shared_topo_vectors.cross.Run_2.x0{Map1,Map2,subj,Hemisphere} = topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_2.Names{Map2}).(Run_2.DTrun).(x_param_2).(variable_2)(:,get_vec_ids_map_2);
                    all_combis.shared_topo_vectors.normals{Map1,Map2,subj,Hemisphere} = topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_1.Names{Map1}).norm_surf(:,get_vec_ids_map_1);
                    if Run_1.Type == "Timing"
                        all_combis.shared_topo_vectors.cross.Run_1.y0{Map1,Map2,subj,Hemisphere} = topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_1.Names{Map1}).(Run_1.DTrun).(y_param_1).(variable_1)(:,get_vec_ids_map_1);
                    end
                    if Run_2.Type == "Timing"
                        all_combis.shared_topo_vectors.cross.Run_2.y0{Map1,Map2,subj,Hemisphere} = topo.(subjNames{subj}).(Hemispheres{Hemisphere}).(Run_2.Names{Map2}).(Run_2.DTrun).(y_param_2).(variable_2)(:,get_vec_ids_map_2);
                    end
                end

                save(all_combi_path,'all_combis')

                do_run_1 = 0;
            end
            do_run_2 = 0;
        end
    end
end

