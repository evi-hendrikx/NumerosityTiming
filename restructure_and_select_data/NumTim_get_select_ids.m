function select_ids = NumTim_get_select_ids(all_combis, whichCombi,DT_runs)
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
%       5 = use voxel selection from other data structure (like 3 but overlap with 4 from another comparison)
% DT_runs: to be compared datatypes
%
% Output
% select_ids: ids indicating which combinations to include in the analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

select_ids = cell(length(DT_runs),1);

%% find indices most overlapping map
amount_overlap = arrayfun(@(x) length(x{:}), all_combis.shared_coord_ids);
% if map doesn't exist
amount_overlap(all_combis.Map1 == "") = NaN;

% remove "All" maps (which automatically have the largest amount of
% overlap)
reduced_amount_overlap = amount_overlap(1:end-1,1:end-1,:,:);


if whichCombi == 2
    % for each DTrun find most overlapping map
    for run = 1:length(DT_runs)
        clear max_overlap id_matrix

        %% compute highest proportions and find index
        % calculate the maximum overlap per map and their indices in the
        % matrix
        % works because you change over which dimension you do it per run
        opposite_dimension = 3 - run;
        [max_overlap, stat.data.(DT_runs{run}).highest_prop_id] = max(reduced_amount_overlap,[],opposite_dimension,'omitnan');

        % NOTE with 0 overlap, first "other map" is chosen. Shouldn't be a problem:
        %     max overlap IS 0, so value max proportion correct
        %     in ANOVA only the main map for which the max proportion is calculated is taken into account (not the other one)
        %     correlation and topography will not be calculated for 0 overlap

        %% transform indices into linear indices
        % make this in the size of all_combi (incl All) again, so it can be
        % indiced in the same way as the other whichCombis
        % this can be done better, but I guess this works
        sz = size(reduced_amount_overlap);
        id_matrix = zeros(size(all_combis.shared_coord_ids));
        if run == 1
            for sz_1 = 1:sz(1) % map
                for sz_3 = 1:sz(3) % subjects
                    for sz_4 = 1:sz(4) % hemispheres
                        if (isnan(max_overlap(sz_1,:,sz_3, sz_4))) % no map
                            continue
                        else
                            id_matrix(sz_1,stat.data.(DT_runs{run}).highest_prop_id(sz_1,:,sz_3, sz_4),sz_3,sz_4) = 1;
                        end
                    end
                end
            end
        else
            for sz_2 = 1:sz(2)
                for sz_3 = 1:sz(3)
                    for sz_4 = 1:sz(4)
                        if (isnan(max_overlap(:,sz_2,sz_3, sz_4))) %no map
                            continue
                        else
                            id_matrix(stat.data.(DT_runs{run}).highest_prop_id(:,sz_2,sz_3, sz_4),sz_2,sz_3,sz_4) = 1;
                        end
                    end
                end
            end
        end
        select_ids{run} = find(id_matrix==1);
    end
else

    %% select appropriate data for whichCombi
    % data per DT_run compared to "All" maps of the other DT_run (1);
    % or the same map (only if both are Timing or both are Numerosity maps) (without all) (3)
    for run = 1:length(DT_runs)

        % Find indices requested combinations
        corresponding_map = strcat('Map',char(num2str(run)));
        other_map = strcat('Map',char(num2str(3-run)));

        if whichCombi == 1 || whichCombi == 4 || whichCombi == 5 % all or all overlapping
            select_ids{run} = contains(all_combis.(other_map),"All") & ~contains(all_combis.(corresponding_map), "All");
        elseif whichCombi == 3  % same
            select_ids{run} = (all_combis.Map1 == all_combis.Map2) & (~contains(all_combis.Map1, "All")) & (all_combis.Map1 ~= "");
        end

    end
end


