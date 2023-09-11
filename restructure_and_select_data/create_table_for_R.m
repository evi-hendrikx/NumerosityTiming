function tbl = create_table_for_R(stat,tbl_path,mean_angles_over_runs)
%% gets topo info and puts it in a suitable format for R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% stat: structure with statistical info
% tbl_path: where you want to store the resulting table
% mean_angles_over_runs (boolean): 1 if mean angle over even (vs odd) and odd (vs even) runs
%
% Output
% tbl: table in correct format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DT_runs = fieldnames(stat.data);
rois = [];
angs = [];

par_1 = fieldnames(stat.data.(DT_runs{1}).topo);
par_2 = fieldnames(stat.data.(DT_runs{1}).topo.(par_1{1}));



for p1 = 1:length(par_1)
    for p2 = 1:length(par_2)
        clear tbl
        rois = [];
        angs = [];

        if mean_angles_over_runs ~= 1
            for run = 1:length(DT_runs)
                mapnames = unique(stat.data.(DT_runs{run}).map);
                for map = 1:length(mapnames)
                    angles = stat.data.(DT_runs{run}).topo.(par_1{p1}).(par_2{p2})(stat.data.(DT_runs{run}).map == mapnames(map))';
                    rois = [rois repmat(mapnames(map),1,length(angles))];
                    angs = [angs angles'];
                end
            end

        else
            mapnames = unique(stat.data.(DT_runs{1}).map);
            for map = 1:length(mapnames)
                ang1 = stat.data.(DT_runs{1}).topo.(par_1{p1}).(par_2{p2})(stat.data.(DT_runs{1}).map == mapnames(map));
                ang2 = stat.data.(DT_runs{2}).topo.(par_1{p1}).(par_2{p2})(stat.data.(DT_runs{1}).map == mapnames(map));

                % make sure means are not NaN if only 1 is NaN
                ang1(isnan(ang1)&~isnan(ang2)) = ang2(isnan(ang1)&~isnan(ang2));
                ang2(isnan(ang2)&~isnan(ang1)) = ang1(isnan(ang2)&~isnan(ang1));

                angles = rad2deg(circ_mean([deg2rad(ang1);deg2rad(ang2)]',[],2));
                rois = [rois repmat(mapnames(map),1,length(angles))];
                angs = [angs angles'];
            end

        end
        tbl = table;
        tbl.map = rois';
        tbl.angles = angs';
        par_tbl_path = insertAfter(tbl_path,'topoTbl', strcat(par_1{p1},'_',par_2{p2},'_meanRuns=',char(num2str(mean_angles_over_runs))));
        writetable(tbl,par_tbl_path)
    end
end
end
