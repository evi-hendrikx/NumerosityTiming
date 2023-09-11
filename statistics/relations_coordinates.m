%% Calculate relations with coordinates (for proportion and topo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% [analysis subject: proportion, topo] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% [steps taken for this subject]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PROP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% general info maps and structure coordinate files
load('center_coords.mat')
load('CenterCoordinatesTiming.mat')
file_general = 'NumerosityAll_TimingAll_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat';
stat_file = ['stat_' file_general];
save_file = ['coords_' file_general];

load(stat_file)

map_types = {'num','time'};


rois.time = ["TLO","TTOP","TTOA","TPO","TLS","TPCI","TPCM","TPCS","TFI","TFS"];
rois.num = ["NLO","NTO","NPO","NPCI","NPCM","NPCS","NFI","NFS"];

num_subj_order_n27 = {'S1','S2','S4','S7','S8','S3','S6','S5'};
num_map_order_n27 = {'','NLO','NLO','NLO','','NTO','NTO','NTO','','NPO','NPO','NPO','','NPCI','NPCI','NPCI','','NPCM','NPCM','NPCM','','NPCS','NPCS','NPCS','','NFI','NFI','NFI','','NFS','NFS','NFS'};
num_subjs = unique(stat.data.NumerosityAll.subj,'stable');
num_hemis = unique(stat.data.NumerosityAll.hemi,'stable');
num_maps = unique(stat.data.NumerosityAll.map,'stable');


time_subj_order_n27 = {'S1','S2','S4','S7','S8','S3','S6','S5'};
time_map_order_n27 = {'','TLO','TLO','TLO','','TTOP','TTOP','TTOP','','TTOA','TTOA','TTOA','','TPO','TPO','TPO','','TLS','TLS','TLS','','TPCI','TPCI','TPCI','','TPCM','TPCM','TPCM','','TPCS','TPCS','TPCS','','TFI','TFI','TFI','','TFS','TFS','TFS'};
time_subjs = unique(stat.data.TimingAll.subj,'stable');
time_hemis = unique(stat.data.TimingAll.hemi,'stable');
time_maps = unique(stat.data.TimingAll.map,'stable');


%% ancovas
pvals_y = [];
pvals_z = [];
pvals_subj = [];
pvals_hemi = [];
for mt = 1:length(map_types)
    if string(map_types{mt}) == "num"
        DT_run = 'NumerosityAll';
    else
        DT_run = 'TimingAll';
    end

    eval(['subjs = ' map_types{mt} '_subjs;'])
    eval(['hemis = ' map_types{mt} '_hemis;'])
    eval(['maps = ' map_types{mt} '_maps;'])
    eval(['map_order_n27 = ' map_types{mt} '_map_order_n27;'])
    eval(['subj_order_n27 = ' map_types{mt} '_subj_order_n27;'])

    % get coords in right order
    y_coords = NaN(size(stat.data.(DT_run).prop)); z_coords = NaN(size(stat.data.(DT_run).prop));
    for subj = 1:length(subjs)
        for hemi = 1:length(hemis)
            for map = 1:length(maps)

                if mt == 1
                    eval(['CenterCoordsN27 = CenterCoords', hemis{hemi}, 'AllN27']);
                    [~,n27_id_map,~] = intersect(map_order_n27, maps{map});
                    [~,n27_id_subj,~] = intersect(subj_order_n27, subjs{subj});

                    y_coords(stat.data.(DT_run).subj == subjs{subj} & stat.data.(DT_run).hemi == hemis{hemi} & stat.data.(DT_run).map == maps{map}) = CenterCoordsN27(n27_id_subj,n27_id_map + 1);
                    z_coords(stat.data.(DT_run).subj == subjs{subj} & stat.data.(DT_run).hemi == hemis{hemi} & stat.data.(DT_run).map == maps{map}) = CenterCoordsN27(n27_id_subj,n27_id_map + 2);
                else
                    eval(['CenterCoordsN27 = CenterCoords', hemis{hemi}, maps{map}]);
                    [~,n27_id_subj,~] = intersect(subj_order_n27, subjs{subj});

                    y_coords(stat.data.(DT_run).subj == subjs{subj} & stat.data.(DT_run).hemi == hemis{hemi} & stat.data.(DT_run).map == maps{map}) = CenterCoordsN27(n27_id_subj,2);
                    z_coords(stat.data.(DT_run).subj == subjs{subj} & stat.data.(DT_run).hemi == hemis{hemi} & stat.data.(DT_run).map == maps{map}) = CenterCoordsN27(n27_id_subj,3);
                end

            end
        end
    end

    covariates = {y_coords,z_coords,stat.data.(DT_run).subj, stat.data.(DT_run).hemi};
    covariate_names = {'y','z','subj','hemi'};
    [~, tbl, stats, ~] = anovan(stat.data.(DT_run).prop, covariates,'Continuous',[1 2],'varnames',covariate_names);
    pvals_subj = [pvals_subj tbl{4,7}];
    pvals_hemi = [pvals_hemi tbl{5,7}];
    pvals_y = [pvals_y tbl{2,7}];
    pvals_z = [pvals_z tbl{3,7}];

    coords_relation.proportion.(DT_run).tbl = tbl;
    coords_relation.proportion.(DT_run).stats = stats;


end

[~,~,~,adj_p_subj] = fdr_bh(pvals_subj);
[~,~,~,adj_p_hemi] = fdr_bh(pvals_hemi);
[~,~,~,adj_p_y] = fdr_bh(pvals_y);
[~,~,~,adj_p_z] = fdr_bh(pvals_z);

for mt = 1:length(map_types)
    if string(map_types{mt}) == "num"
        DT_run = 'NumerosityAll';
    else
        DT_run = 'TimingAll';
    end
    coords_relation.proportion.(DT_run).adj_p.subj =adj_p_subj(mt);
    coords_relation.proportion.(DT_run).adj_p.hemi =adj_p_hemi(mt);
    coords_relation.proportion.(DT_run).adj_p.y =adj_p_y(mt);
    coords_relation.proportion.(DT_run).adj_p.z =adj_p_z(mt);
end

%% TOPO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% general info
clearvars -except coords_relation
close all

load('center_coords.mat')
load('CenterCoordinatesTiming.mat')

file_general = 'NumerosityAll_TimingAll_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat';
stat_file = ['stat_' file_general];
save_file = ['coords_' file_general];

load(stat_file)

map_types = {'num','time'};

%% make mean Timing map coordinates (topo stats are on summary statistics over participants)
rois.time = ["TLO","TTOP","TTOA","TPO","TLS","TPCI","TPCM","TPCS","TFI","TFS"];
for roi = 1:length(rois.time)
    eval(['mean_coords = nanmean([CenterCoordsLeft', char(rois.time(roi)), '; CenterCoordsRight', char(rois.time(roi)), ']);'])
    y_coords.time(roi) = mean_coords(2);
    z_coords.time(roi) = mean_coords(3);
end

% make mean Numerosity map coordinates
rois.num = ["NLO","NTO","NPO","NPCI","NPCM","NPCS","NFI","NFS"];
num_map_order_n27 = {'','NLO','NLO','NLO','','NTO','NTO','NTO','','NPO','NPO','NPO','','NPCI','NPCI','NPCI','','NPCM','NPCM','NPCM','','NPCS','NPCS','NPCS','','NFI','NFI','NFI','','NFS','NFS','NFS'};
for roi = 1:length(rois.num)
    [~,n27_id_map,~] = intersect(num_map_order_n27, rois.num(roi));

    y_coords.num(roi) = nanmean([CenterCoordsLeftAllN27(:,n27_id_map + 1); CenterCoordsRightAllN27(:,n27_id_map + 1)]);
    z_coords.num(roi) = nanmean([CenterCoordsLeftAllN27(:,n27_id_map + 2); CenterCoordsRightAllN27(:,n27_id_map + 2)]);
end

%% bayes factors (taken from previous analyes
BF.NumerosityAll.x0.x0 = [0.3358220,3.7892062,0.8212109,0.3130050,0.4090667,4.2467176,1.6032793, 1.6517540];
BF.NumerosityAll.x0.y0 = [0.4700467,0.6256568,0.7962206,0.2606131,0.4692343,282.1445297,0.9649128,2.2703932];
BF.TimingAll.x0.x0 = [0.3358220,3.7892062,0.6948372,0.2707176,0.2695510,1.4263719,1.6032793,0.3305964,0.6731898];
BF.TimingAll.x0.y0 = [0.4700467,0.6256568,0.6788992,0.2612219,0.2614731,44.2828100,0.9649128,0.3422824,1.8632803];
BF_rois.NumerosityAll = {'NFI','NFS','NLO','NPCI','NPCM','NPCS','NPO','NTO'};
BF_rois.TimingAll = {'TFI','TFS','TLO','TPCI','TPCM','TPCS','TPO','TTOA','TTOP'};

%% partial correlations
pvals_range_1 = [];
pvals_range_2 = [];
for mt = 1:length(map_types)
    if string(map_types{mt}) == "num"
        DT_run = 'NumerosityAll';
    else
        DT_run = 'TimingAll';
    end
    [roi_order,new_ids_rois,new_ids_bf_order]=intersect(rois.(map_types{mt}),BF_rois.(DT_run),'stable');

    [r,pvals] = partialcorr([BF.(DT_run).x0.x0(new_ids_bf_order)' y_coords.(map_types{mt})(new_ids_rois)' z_coords.(map_types{mt})(new_ids_rois)'],'type','Spearman','rows','complete');
    coords_relation.topo.(DT_run).x0x0.y.p = pvals(1,2);
    coords_relation.topo.(DT_run).x0x0.z.p = pvals(1,3);
    coords_relation.topo.(DT_run).x0x0.y.r = r(1,2);
    coords_relation.topo.(DT_run).x0x0.z.r = r(1,3);

    pvals_range_1 = [pvals_range_1 coords_relation.topo.(DT_run).x0x0.y.p];
    pvals_range_2 = [pvals_range_2 coords_relation.topo.(DT_run).x0x0.z.p];

    [r,pvals] = partialcorr([BF.(DT_run).x0.y0(new_ids_bf_order)' y_coords.(map_types{mt})(new_ids_rois)' z_coords.(map_types{mt})(new_ids_rois)'],'type','Spearman','rows','complete');
    coords_relation.topo.(DT_run).x0y0.y.p = pvals(1,2);
    coords_relation.topo.(DT_run).x0y0.z.p = pvals(1,3);
    coords_relation.topo.(DT_run).x0y0.y.r = r(1,2);
    coords_relation.topo.(DT_run).x0y0.z.r = r(1,3);

    pvals_range_1 = [pvals_range_1 coords_relation.topo.(DT_run).x0y0.y.p];
    pvals_range_2 = [pvals_range_2 coords_relation.topo.(DT_run).x0y0.z.p];


end

[~,~,~,adj_p_range_1] = fdr_bh(pvals_range_1);
[~,~,~,adj_p_range_2] = fdr_bh(pvals_range_2);

counter = 1;
for mt = 1:length(map_types)
    if string(map_types{mt}) == "num"
        DT_run = 'NumerosityAll';
    else
        DT_run = 'TimingAll';
    end

    coords_relation.topo.(DT_run).x0x0.y.adj_p = adj_p_range_1(counter);
    coords_relation.topo.(DT_run).x0x0.z.adj_p = adj_p_range_2(counter);

    counter = counter + 1;

    coords_relation.topo.(DT_run).x0y0.y.adj_p = adj_p_range_1(counter);
    coords_relation.topo.(DT_run).x0y0.z.adj_p = adj_p_range_2(counter);

    counter = counter + 1;
end

save(save_file, 'coords_relation')