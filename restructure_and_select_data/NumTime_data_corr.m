function stat = NumTime_data_corr(stat,whichCombi,all_combis,select_ids,DT_runs,run,minAmountCorrelation, selection_all_combis,selection_DT_runs,selection_run)
%% Selects and computes correlations from requested ids and assigns them to stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% stat: structure containing statistical info from all the performed tests
% whichCombi: 1 = compare to "All" maps; 2 = compare to most overlapping
%       map; 3 = compare to same map; 4 = compare to all overlapping maps
%       5 = use voxel selection from other data structure
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
% minAmountCorrelation: minimum amount of voxels a map should have in order
%       to calculate the correlation coefficients (and topo angles)
%
% Output
% stat: structure containing statistical info from all the performed tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get preferred values of shared voxels
if whichCombi ~=4 && whichCombi ~= 5
    values_1_x = all_combis.preferred_values.Run_1.x0(select_ids);
    values_2_x = all_combis.preferred_values.Run_2.x0(select_ids);
    if contains(DT_runs{1}, 'Timing')
        values_1_y = all_combis.preferred_values.Run_1.y0(select_ids);
    end
    if contains(DT_runs{2}, 'Timing')
        values_2_y = all_combis.preferred_values.Run_2.y0(select_ids);
    end

else

    for datapoint = 1:length(stat.data.(DT_runs{run}).subj)
        % How sure am I this data points is exactly the same order?
        % pretty sure: difference is kicking out all maps. But all of these
        % should only yield one value per data point, so that order cannot
        % get mixed up
        eval(['values_1_x{datapoint} = [all_combis.preferred_values.Run_1.x0{all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & all_combis.Map',char(num2str(run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(all_combis.Map',char(num2str(3-run)),', "All")}];'])
        eval(['values_2_x{datapoint} = [all_combis.preferred_values.Run_2.x0{all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & all_combis.Map',char(num2str(run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(all_combis.Map',char(num2str(3-run)),', "All")}];'])
        if contains(DT_runs{1}, 'Timing')
            eval(['values_1_y{datapoint} = [all_combis.preferred_values.Run_1.y0{all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & all_combis.Map',char(num2str(run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(all_combis.Map',char(num2str(3-run)),', "All")}];'])
        end
        if contains(DT_runs{2}, 'Timing')
            eval(['values_2_y{datapoint} = [all_combis.preferred_values.Run_2.y0{all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & all_combis.Map',char(num2str(run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(all_combis.Map',char(num2str(3-run)),', "All")}];'])
        end
    end

    if whichCombi == 5
        % also get values of the voxels used for selection
        for datapoint = 1:length(stat.data.(DT_runs{run}).subj)
            eval(['shared_coords{datapoint} = [all_combis.shared_coord_ids{all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & all_combis.Map',char(num2str(run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(all_combis.Map',char(num2str(3-run)),', "All")}];'])

            eval(['selection_shared_coords{datapoint} = [selection_all_combis.shared_coord_ids{selection_all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & selection_all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & selection_all_combis.Map',char(num2str(selection_run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(selection_all_combis.Map',char(num2str(3-selection_run)),', "All")}];'])
            eval(['selection_values_1_x{datapoint}= [selection_all_combis.preferred_values.Run_1.x0{selection_all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & selection_all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & selection_all_combis.Map',char(num2str(selection_run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(selection_all_combis.Map',char(num2str(3-selection_run)),', "All")}];'])
            eval(['selection_values_2_x{datapoint} = [selection_all_combis.preferred_values.Run_2.x0{selection_all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & selection_all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & selection_all_combis.Map',char(num2str(selection_run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(selection_all_combis.Map',char(num2str(3-selection_run)),', "All")}];'])
            if  contains(selection_DT_runs{1},'Timing')
                eval(['selection_values_1_y{datapoint} = [selection_all_combis.preferred_values.Run_1.y0{selection_all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & selection_all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & selection_all_combis.Map',char(num2str(selection_run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(selection_all_combis.Map',char(num2str(3-selection_run)),', "All")}];'])
            end
            if  contains(selection_DT_runs{2},'Timing')
                eval(['selection_values_2_y{datapoint} = [selection_all_combis.preferred_values.Run_2.y0{selection_all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & selection_all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & selection_all_combis.Map',char(num2str(selection_run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(selection_all_combis.Map',char(num2str(3-selection_run)),', "All")}];'])
            end
        end

        for  datapoint = 1:length(stat.data.(DT_runs{run}).subj)

            % within the area that was overlapping for the other selection
            % Which voxels also overlap for these two runs?
            [~,id_overlap_voxels{datapoint},selection_id_overlap_voxels{datapoint}] = intersectCols(shared_coords{datapoint},selection_shared_coords{datapoint});

            % voxels within overlapping area shared by all
            values_1_x{datapoint} = values_1_x{datapoint}(id_overlap_voxels{datapoint});
            if contains(DT_runs{1}, 'Timing')
                values_1_y{datapoint} = values_1_y{datapoint}(id_overlap_voxels{datapoint});
            end
            values_2_x{datapoint} = values_2_x{datapoint}(id_overlap_voxels{datapoint});
            if contains(DT_runs{2}, 'Timing')
                values_2_y{datapoint} = values_2_y{datapoint}(id_overlap_voxels{datapoint});
            end

            selection_values_1_x{datapoint} = selection_values_1_x{datapoint}(selection_id_overlap_voxels{datapoint});
            selection_values_2_x{datapoint} = selection_values_2_x{datapoint}(selection_id_overlap_voxels{datapoint});
            if contains(selection_DT_runs{1}, 'Timing')
                selection_values_1_y{datapoint} = selection_values_1_y{datapoint}(selection_id_overlap_voxels{datapoint});
            end
            if contains(selection_DT_runs{2}, 'Timing')
                selection_values_2_y{datapoint} = selection_values_2_y{datapoint}(selection_id_overlap_voxels{datapoint});
            end
        end
    end

end

%% calulate correlations
% we want to know whether there is a relationship (not necessarily a linear one), so Spearman correlation


for datapoint = 1:length(stat.data.(DT_runs{run}).subj)
    if length(values_1_x{datapoint}) >= minAmountCorrelation

    

        % Numerosity only has preferred numerosity value (x0) and timing
        % has preferred duration (x0) and preferred period (y0)
        stat.data.(DT_runs{run}).correlation.x0.x0(datapoint) = corr(values_1_x{datapoint}', values_2_x{datapoint}', 'Type', 'Spearman','rows','pairwise'); % pairwise shouldn't make a difference here as we only put in two colums

            stat.data.(DT_runs{run}).preferred_values.std.Run_1_x0(datapoint) = std(values_1_x{datapoint}', 'omitnan');
            stat.data.(DT_runs{run}).preferred_values.iqr.Run_1_x0(datapoint) = iqr(values_1_x{datapoint}');
            stat.data.(DT_runs{run}).preferred_values.std.Run_2_x0(datapoint) = std(values_2_x{datapoint}', 'omitnan');
            stat.data.(DT_runs{run}).preferred_values.iqr.Run_2_x0(datapoint) = iqr(values_2_x{datapoint}');

        if contains(DT_runs{1}, 'Timing')
            stat.data.(DT_runs{run}).correlation.y0.x0(datapoint) = corr(values_1_y{datapoint}', values_2_x{datapoint}', 'Type', 'Spearman','rows','pairwise');
            stat.data.(DT_runs{run}).preferred_values.std.Run_1_y0(datapoint) = std(values_1_y{datapoint}', 'omitnan');
            stat.data.(DT_runs{run}).preferred_values.iqr.Run_1_y0(datapoint) = iqr(values_1_y{datapoint}');
        end
        if contains(DT_runs{2}, 'Timing')
            stat.data.(DT_runs{run}).correlation.x0.y0(datapoint) = corr(values_1_x{datapoint}', values_2_y{datapoint}', 'Type', 'Spearman','rows','pairwise');
            stat.data.(DT_runs{run}).preferred_values.std.Run_2_y0(datapoint) = std(values_2_y{datapoint}', 'omitnan');
            stat.data.(DT_runs{run}).preferred_values.iqr.Run_2_y0(datapoint) = iqr(values_2_y{datapoint}');
        end
        if contains(DT_runs{1}, 'Timing') && contains(DT_runs{2}, 'Timing')
            stat.data.(DT_runs{run}).correlation.y0.y0(datapoint) = corr(values_1_y{datapoint}', values_2_y{datapoint}', 'Type', 'Spearman','rows','pairwise');
        end
    else
        stat.data.(DT_runs{run}).correlation.x0.x0(datapoint) = NaN;
        stat.data.(DT_runs{run}).preferred_values.std.Run_1_x0(datapoint) = NaN;
        stat.data.(DT_runs{run}).preferred_values.std.Run_2_x0(datapoint) = NaN;
        stat.data.(DT_runs{run}).preferred_values.iqr.Run_1_x0(datapoint) = NaN;
        stat.data.(DT_runs{run}).preferred_values.iqr.Run_2_x0(datapoint) = NaN;
        if contains(DT_runs{1}, 'Timing')
            stat.data.(DT_runs{run}).correlation.y0.x0(datapoint) = NaN;
            stat.data.(DT_runs{run}).preferred_values.std.Run_1_y0(datapoint) = NaN;
            stat.data.(DT_runs{run}).preferred_values.iqr.Run_1_y0(datapoint) = NaN;
        end
        if contains(DT_runs{2}, 'Timing')
            stat.data.(DT_runs{run}).correlation.x0.y0(datapoint) = NaN;
            stat.data.(DT_runs{run}).preferred_values.std.Run_2_y0(datapoint) = NaN;
            stat.data.(DT_runs{run}).preferred_values.iqr.Run_2_y0(datapoint) = NaN;
        end
        if contains(DT_runs{1}, 'Timing') && contains(DT_runs{2}, 'Timing')
            stat.data.(DT_runs{run}).correlation.y0.y0(datapoint) = NaN;
        end
    end


    if whichCombi == 5
        % also calculate correlations for other voxel selection
        if length(selection_values_1_x{datapoint}) >= minAmountCorrelation

            % Numerosity only has preferred numerosity value (x0) and timing
            % has preferred duration (x0) and preferred period (y0)
            stat.data.(selection_DT_runs{run}).selection_correlation.x0.x0(datapoint) = corr(selection_values_1_x{datapoint}', selection_values_2_x{datapoint}', 'Type', 'Spearman','rows','pairwise'); % pairwise shouldn't make a difference here as we only put in two colums

            if contains(selection_DT_runs{1}, 'Timing')
                stat.data.(selection_DT_runs{run}).selection_correlation.y0.x0(datapoint) = corr(selection_values_1_y{datapoint}', selection_values_2_x{datapoint}', 'Type', 'Spearman','rows','pairwise');
            end
            if contains(selection_DT_runs{2}, 'Timing')
                stat.data.(selection_DT_runs{run}).selection_correlation.x0.y0(datapoint) = corr(selection_values_1_x{datapoint}', selection_values_2_y{datapoint}', 'Type', 'Spearman','rows','pairwise');
            end
            if contains(selection_DT_runs{1}, 'Timing') && contains(selection_DT_runs{2}, 'Timing')
                stat.data.(selection_DT_runs{run}).selection_correlation.y0.y0(datapoint) = corr(selection_values_1_y{datapoint}', selection_values_2_y{datapoint}', 'Type', 'Spearman','rows','pairwise');
            end

        else

            stat.data.(selection_DT_runs{run}).selection_correlation.x0.x0(datapoint) = NaN;
            if contains(selection_DT_runs{1}, 'Timing')
                stat.data.(selection_DT_runs{run}).selection_correlation.y0.x0(datapoint) = NaN;
            end
            if contains(selection_DT_runs{2}, 'Timing')
                stat.data.(selection_DT_runs{run}).selection_correlation.x0.y0(datapoint) = NaN;
            end
            if contains(selection_DT_runs{1}, 'Timing') && contains(selection_DT_runs{2}, 'Timing')
                stat.data.(selection_DT_runs{run}).selection_correlation.y0.y0(datapoint) = NaN;
            end

        end
    end
end


