%% Calculate relations with range of values within the maps (for correlation and topo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% [analysis subject: correlation, topo] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% [steps taken for this subject]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CORR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% general info

file_general = 'NumerosityAll_TimingAll_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat';
stat_file = ['stat_' file_general];
save_file = ['corr_range_' file_general];

load(stat_file)

map_types = {'num','time'};

%% rois
rois.time = ["TLO","TTOP","TTOA","TPO","TLS","TPCI","TPCM","TPCS","TFI","TFS"];
rois.num = ["NLO","NTO","NPO","NPCI","NPCM","NPCS","NFI","NFS"];

%% partial correlations

pvals_range_1 = [];
pvals_range_2 = [];
for mt = 1:length(map_types)
    if string(map_types{mt}) == "num"
        DT_run = 'NumerosityAll';
    else
        DT_run = 'TimingAll';
    end

    [r,pvals] = partialcorr([abs(stat.data.(DT_run).correlation.x0.x0)' stat.data.(DT_run).preferred_values.iqr.Run_1_x0' stat.data.(DT_run).preferred_values.iqr.Run_2_x0'],'type','Spearman','rows','complete');
    par_corr.(DT_run).x0x0.var_1.p = pvals(1,2);
    par_corr.(DT_run).x0x0.var_2.p = pvals(1,3);
    par_corr.(DT_run).x0x0.var_1.r = r(1,2);
    par_corr.(DT_run).x0x0.var_2.r = r(1,3);

    pvals_range_1 = [pvals_range_1 par_corr.(DT_run).x0x0.var_1.p];
    pvals_range_2 = [pvals_range_2 par_corr.(DT_run).x0x0.var_2.p];

    [r,pvals] = partialcorr([abs(stat.data.(DT_run).correlation.x0.y0)' stat.data.(DT_run).preferred_values.iqr.Run_1_x0' stat.data.(DT_run).preferred_values.iqr.Run_2_y0'],'type','Spearman','rows','complete');
    par_corr.(DT_run).x0y0.var_1.p = pvals(1,2);
    par_corr.(DT_run).x0y0.var_2.p = pvals(1,3);
    par_corr.(DT_run).x0y0.var_1.r = r(1,2);
    par_corr.(DT_run).x0y0.var_2.r = r(1,3);

    pvals_range_1 = [pvals_range_1 par_corr.(DT_run).x0y0.var_1.p];
    pvals_range_2 = [pvals_range_2 par_corr.(DT_run).x0y0.var_2.p];
end

[~,~,~,adj_p_range_1] = fdr_bh(pvals_range_1);
[~,~,~,adj_p_range_2] = fdr_bh(pvals_range_2);

counter = 0;
for mt = 1:length(map_types)

    if string(map_types{mt}) == "num"
        DT_run = 'NumerosityAll';
    else
        DT_run = 'TimingAll';
    end

    counter = counter + 1;

    par_corr.(DT_run).x0x0.var_1.adj_p = adj_p_range_1(counter);
    par_corr.(DT_run).x0x0.var_2.adj_p = adj_p_range_2(counter);

    counter = counter + 1;

    par_corr.(DT_run).x0y0.var_1.adj_p = adj_p_range_1(counter);
    par_corr.(DT_run).x0y0.var_2.adj_p = adj_p_range_2(counter);


end

save(save_file,'par_corr')

%% TOPO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except par_corr
close all

%% general info
file_general = 'NumerosityAll_TimingAll_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat';
stat_file = ['stat_' file_general];
save_file = ['corr_range_' file_general];

load(stat_file)

map_types = {'num','time'};

rois.time = ["TLO","TTOP","TTOA","TPO","TLS","TPCI","TPCM","TPCS","TFI","TFS"];
rois.num = ["NLO","NTO","NPO","NPCI","NPCM","NPCS","NFI","NFS"];


%% BFs (copied from previous analyses)
BF.NumerosityAll.x0.x0 = [0.3358220,3.7892062,0.8212109,0.3130050,0.4090667,4.2467176,1.6032793, 1.6517540];
BF.NumerosityAll.x0.y0 = [0.4700467,0.6256568,0.7962206,0.2606131,0.4692343,282.1445297,0.9649128,2.2703932];
BF.TimingAll.x0.x0 = [0.3358220,3.7892062,0.6948372,0.2707176,0.2695510,1.4263719,1.6032793,0.3305964,0.6731898];
BF.TimingAll.x0.y0 = [0.4700467,0.6256568,0.6788992,0.2612219,0.2614731,44.2828100,0.9649128,0.3422824,1.8632803];
BF_rois.NumerosityAll = {'NFI','NFS','NLO','NPCI','NPCM','NPCS','NPO','NTO'};
BF_rois.TimingAll = {'TFI','TFS','TLO','TPCI','TPCM','TPCS','TPO','TTOA','TTOP'};


%% partial correlations with iqr preferred values of the maps
pvals_range_1 = [];
pvals_range_2 = [];

for mt = 1:length(map_types)
    if string(map_types{mt}) == "num"
        DT_run = 'NumerosityAll';
    else
        DT_run = 'TimingAll';
    end

    mean_iqr_1_x0 = [];
    mean_iqr_2_x0 = [];
    mean_iqr_2_y0 = [];

    for roi = 1:length(BF_rois.(DT_run))
        try
            [~,BF_topo.(DT_run).x0x0.range_1.pnorm.(BF_rois.(DT_run){roi})] = force_swtest(stat.data.(DT_run).preferred_values.iqr.Run_1_x0(string(stat.data.(DT_run).map) == BF_rois.(DT_run){roi}));
        catch
            BF_topo.(DT_run).x0x0.range_1.pnorm.(BF_rois.(DT_run){roi}) = NaN;
        end
        try
            [~,BF_topo.(DT_run).x0x0.range_2.pnorm.(BF_rois.(DT_run){roi})] = force_swtest(stat.data.(DT_run).preferred_values.iqr.Run_2_x0(string(stat.data.(DT_run).map) == BF_rois.(DT_run){roi}));
        catch
            BF_topo.(DT_run).x0x0.range_2.pnorm.(BF_rois.(DT_run){roi}) = NaN;
        end
        try
            [~,BF_topo.(DT_run).x0y0.range_2.pnorm.(BF_rois.(DT_run){roi})] = force_swtest(stat.data.(DT_run).preferred_values.iqr.Run_2_y0(string(stat.data.(DT_run).map) == BF_rois.(DT_run){roi}));
        catch
            BF_topo.(DT_run).x0y0.range_2.pnorm.(BF_rois.(DT_run){roi}) = NaN;
        end

        mean_iqr_1_x0 = [mean_iqr_1_x0 nanmedian(stat.data.(DT_run).preferred_values.iqr.Run_1_x0(string(stat.data.(DT_run).map) == BF_rois.(DT_run){roi}))];
        mean_iqr_2_x0 = [mean_iqr_2_x0 nanmedian(stat.data.(DT_run).preferred_values.iqr.Run_2_x0(string(stat.data.(DT_run).map) == BF_rois.(DT_run){roi}))];
        mean_iqr_2_y0 = [mean_iqr_2_y0 nanmedian(stat.data.(DT_run).preferred_values.iqr.Run_2_y0(string(stat.data.(DT_run).map) == BF_rois.(DT_run){roi}))];
    end

    [r,p] = partialcorr([BF.(DT_run).x0.x0' mean_iqr_1_x0' mean_iqr_2_x0'],'type','Spearman','rows','complete');
    BF_topo.(DT_run).x0x0.range_1.r = r(1,2);
    BF_topo.(DT_run).x0x0.range_1.p = p(1,2);
    BF_topo.(DT_run).x0x0.range_2.r = r(1,3);
    BF_topo.(DT_run).x0x0.range_2.p = p(1,3);

    pvals_range_1 = [pvals_range_1 BF_topo.(DT_run).x0x0.range_1.p];
    pvals_range_2 = [pvals_range_2 BF_topo.(DT_run).x0x0.range_2.p];


    [r,p] = partialcorr([BF.(DT_run).x0.y0' mean_iqr_1_x0' mean_iqr_2_y0'],'type','Spearman','rows','complete');
    BF_topo.(DT_run).x0y0.range_1.r = r(1,2);
    BF_topo.(DT_run).x0y0.range_1.p = p(1,2);
    BF_topo.(DT_run).x0y0.range_2.r = r(1,3);
    BF_topo.(DT_run).x0y0.range_2.p = p(1,3);
    pvals_range_1 = [pvals_range_1 BF_topo.(DT_run).x0y0.range_1.p];
    pvals_range_2 = [pvals_range_2 BF_topo.(DT_run).x0y0.range_2.p];
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

    BF_topo.(DT_run).x0x0.range_1.adj_p = adj_p_range_1(counter);
    BF_topo.(DT_run).x0x0.range_2.adj_p = adj_p_range_2(counter);

    counter = counter + 1;

    BF_topo.(DT_run).x0y0.range_1.adj_p = adj_p_range_1(counter);
    BF_topo.(DT_run).x0y0.range_2.adj_p = adj_p_range_2(counter);

    counter = counter + 1;

end

save(save_file,'par_corr','range_topo','BF_topo')