function stat = NumTime_data_topo(stat,whichCombi,all_combis,select_ids,DT_runs,run,minAmountCorrelation,topo_measurement_per_map,type_angle,stat_path)
%% Selects and computes correlations from requested ids and assigns them to stats
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
% minAmountCorrelation: minimum amount of voxels a map should have in order
%       to calculate the correlation coefficients (and topo angles)
% topo_measurement_per_map: calculate mean or median angle between vectors
% type_angle: angles between topographic vectors can be computed taking
%       into account that the hemispheres are mirror images of each other
%       (medial-lateral), or using the exact same method for both
%       hemispheres (left-right)
% stat_path: path where you want stat (and topo angle table) to be saved
%
% Output
% stat: structure containing statistical info from all the performed tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if whichCombi == 1
    disp("No topo done for All map")
    return
end

%% get values of shared voxels
if whichCombi ~=4
    all_normals = all_combis.shared_topo_vectors.normals(select_ids);
    cross_1_x = all_combis.shared_topo_vectors.cross.Run_1.x0(select_ids);
    cross_2_x = all_combis.shared_topo_vectors.cross.Run_2.x0(select_ids);
    if contains(DT_runs{1}, 'Timing')
        cross_1_y = all_combis.shared_topo_vectors.cross.Run_1.y0(select_ids);
    end
    if contains(DT_runs{2}, 'Timing')
        cross_2_y = all_combis.shared_topo_vectors.cross.Run_2.y0(select_ids);
    end

else
    for datapoint = 1:length(stat.data.(DT_runs{run}).subj)
        eval(['all_normals{datapoint} = [all_combis.shared_topo_vectors.normals{all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & all_combis.Map',char(num2str(run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(all_combis.Map',char(num2str(3-run)),', "All")}];'])
        eval(['cross_1_x{datapoint} = [all_combis.shared_topo_vectors.cross.Run_1.x0{all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & all_combis.Map',char(num2str(run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(all_combis.Map',char(num2str(3-run)),', "All")}];'])
        eval(['cross_2_x{datapoint} = [all_combis.shared_topo_vectors.cross.Run_2.x0{all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & all_combis.Map',char(num2str(run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(all_combis.Map',char(num2str(3-run)),', "All")}];'])
        if contains(DT_runs{1}, 'Timing')
            eval(['cross_1_y{datapoint} = [all_combis.shared_topo_vectors.cross.Run_1.y0{all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & all_combis.Map',char(num2str(run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(all_combis.Map',char(num2str(3-run)),', "All")}];'])
        end
        if contains(DT_runs{2}, 'Timing')
            eval(['cross_2_y{datapoint} = [all_combis.shared_topo_vectors.cross.Run_2.y0{all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & all_combis.Map',char(num2str(run)),' == stat.data.(DT_runs{run}).map(datapoint) & ~contains(all_combis.Map',char(num2str(3-run)),', "All")}];'])
        end
    end
end

%% compute mean vector angle
% select summary statistic to use for each map
if string(topo_measurement_per_map) == "median"
    get_summary_stat = 'circ_median';
    test_input = '2';
elseif string(topo_measurement_per_map) == "mean"
    get_summary_stat = 'circ_mean';
    test_input = '[],2';
end

for datapoint = 1:length(stat.data.(DT_runs{run}).subj)

    %% calulate angles between vectors per voxel
    if length(cross_1_x{datapoint}) >= minAmountCorrelation

        x0x0 = NaN(1,length(cross_1_x{datapoint}));
        normals = all_normals{datapoint};
        cross.x0_1 = cross_1_x{datapoint};
        cross.x0_2 = cross_2_x{datapoint};
        if contains(DT_runs{1}, 'Timing')
            y0x0 = x0x0;
            cross.y0_1 = cross_1_y{datapoint};
        end
        if contains(DT_runs{2}, 'Timing')
            x0y0 = x0x0;
            cross.y0_2 = cross_2_y{datapoint};
        end
        if contains(DT_runs{1}, 'Timing') && contains(DT_runs{2}, 'Timing')
            y0y0 = x0x0;
        end

        % I run vecangle360 per voxel, otherwise it does not do
        % element-wise operations
        for shared_vox = 1:length(cross_1_x{datapoint})
            try
                % you could expect the angles in the different hemispheres to be mirror images
                % of each other
                % angle anticlockwise away from run 1
                if (string(type_angle) == "medial-lateral" && stat.data.(DT_runs{run}).hemi(datapoint) == "Right") || string(type_angle) == "left-right"
                    x0x0(shared_vox) = vecangle360(cross.x0_1(:,shared_vox),cross.x0_2(:,shared_vox),normals(:,shared_vox));
                    if contains(DT_runs{1}, 'Timing')
                        y0x0(shared_vox) = vecangle360(cross.y0_1(:,shared_vox),cross.x0_2(:,shared_vox),normals(:,shared_vox));
                    end
                    if contains(DT_runs{2}, 'Timing')
                        x0y0(shared_vox) = vecangle360(cross.x0_1(:,shared_vox),cross.y0_2(:,shared_vox),normals(:,shared_vox));
                    end
                    if contains(DT_runs{1}, 'Timing') && contains(DT_runs{2}, 'Timing')
                        y0y0(shared_vox) = vecangle360(cross.y0_1(:,shared_vox),cross.y0_2(:,shared_vox),normals(:,shared_vox));
                    end
                    % angle anticlockwise away from run 2
                elseif string(type_angle) == "medial-lateral" && stat.data.(DT_runs{run}).hemi(datapoint) == "Left"
                    x0x0(shared_vox) = vecangle360(cross.x0_2(:,shared_vox),cross.x0_1(:,shared_vox),normals(:,shared_vox));
                    if contains(DT_runs{1}, 'Timing')
                        y0x0(shared_vox) = vecangle360(cross.x0_2(:,shared_vox),cross.y0_1(:,shared_vox),normals(:,shared_vox));
                    end
                    if contains(DT_runs{2}, 'Timing')
                        x0y0(shared_vox) = vecangle360(cross.y0_2(:,shared_vox),cross.x0_1(:,shared_vox),normals(:,shared_vox));
                    end
                    if contains(DT_runs{1}, 'Timing') && contains(DT_runs{2}, 'Timing')
                        y0y0(shared_vox) = vecangle360(cross.y0_2(:,shared_vox),cross.y0_1(:,shared_vox),normals(:,shared_vox));
                    end
                end
            catch
                disp(stat.data.(DT_runs{run}).subj(datapoint))
                disp(stat.data.(DT_runs{run}).hemi(datapoint))
                disp(stat.data.(DT_runs{run}).map(datapoint))
                error("vector has NaN values. Something may have gone wrong in computing shared voxels (these should be above the threshold for both)")
            end

        end

        %% calculate mean angle
        %% and STDEV
        % "Interestingly, multiple quantities have been introduced as analogues to the linear standard deviation. First, the angular deviation is defined ass = sqrt(2 * (1 − R)). This quantity lies in the interval [0, sqrt(2)].
        % Alternatively, the circular standard deviation is defined as s0 = sqrt(−2 * ln(R)) and ranges from 0 to ∞.
        % Generally, the first measure is preferred, as it is bounded, but the two measures deviate little (Zar (1999))."

        stat.data.(DT_runs{run}).topo.x0.x0(datapoint) = eval(['rad2deg(',get_summary_stat, '(deg2rad(x0x0),',test_input,'))']);
        stat.data.(DT_runs{run}).topo_std.x0.x0(datapoint) = rad2deg(circ_std(deg2rad(x0x0)));
        if contains(DT_runs{1}, 'Timing')
            stat.data.(DT_runs{run}).topo.y0.x0(datapoint) = eval(['rad2deg(',get_summary_stat, '(deg2rad(y0x0),',test_input,'))']);
            stat.data.(DT_runs{run}).topo_std.y0.x0(datapoint) = rad2deg(circ_std(deg2rad(y0x0)));
        end
        if contains(DT_runs{2}, 'Timing')
            stat.data.(DT_runs{run}).topo.x0.y0(datapoint) = eval(['rad2deg(',get_summary_stat, '(deg2rad(x0y0),',test_input,'))']);
            stat.data.(DT_runs{run}).topo_std.x0.y0(datapoint) = rad2deg(circ_std(deg2rad(x0y0)));
        end
        if contains(DT_runs{1}, 'Timing') && contains(DT_runs{2}, 'Timing')
            stat.data.(DT_runs{run}).topo.y0.y0(datapoint) = eval(['rad2deg(',get_summary_stat, '(deg2rad(y0y0),',test_input,'))']);
            stat.data.(DT_runs{run}).topo_std.y0.y0(datapoint) = rad2deg(circ_std(deg2rad(y0y0)));
        end
    else
        stat.data.(DT_runs{run}).topo.x0.x0(datapoint) = NaN;
        stat.data.(DT_runs{run}).topo_std.x0.x0(datapoint) = NaN;
        if contains(DT_runs{1}, 'Timing')
            stat.data.(DT_runs{run}).topo.y0.x0(datapoint) = NaN;
            stat.data.(DT_runs{run}).topo_std.y0.x0(datapoint) = NaN;
        end
        if contains(DT_runs{2}, 'Timing')
            stat.data.(DT_runs{run}).topo.x0.y0(datapoint) = NaN;
            stat.data.(DT_runs{run}).topo_std.x0.y0(datapoint) = NaN;
        end
        if contains(DT_runs{1}, 'Timing') && contains(DT_runs{2}, 'Timing')
            stat.data.(DT_runs{run}).topo.y0.y0(datapoint) = NaN;
            stat.data.(DT_runs{run}).topo_std.y0.y0(datapoint) = NaN;
        end
    end

end

end

