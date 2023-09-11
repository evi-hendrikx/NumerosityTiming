%% gets statistical values from matlab structs and outputs them as a table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% different parts of this script were run separately
% %% [analysis subject: proportion, correlation, topo] %%%%%%%%%%%%%%%%%%%%
% %% [different analyses for this subject]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% proportion overlap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
whichCombi = '4';
file_general = ['NumerosityAll_TimingAll_minVEtime=0.2_minVEnum=0.3_whichCombi=' whichCombi '_topo_measure=mean'];
stat_file = ['stat_' file_general '.mat'];
general_save_file = ['data_tbl_' file_general '.csv'];

load(stat_file)
% load plot values (including indicated means and medians)
load(['minVEtime=0.2_minVEnum=0.3_whichCombi=',whichCombi,'_topo_NumerosityAll_TimingAll_prop=1_corr=1_topo=1_plot_values.mat'])


map_types = {'num','time'};

%% rois
rois.time = ["TLO","TTOP","TTOA","TPO","TLS","TPCI","TPCM","TPCS","TFI","TFS"];
rois.num = ["NLO","NTO","NPO","NPCI","NPCM","NPCS","NFI","NFS"];


%% different from chance
% Different from chance is on medians while
% ANOVAs are on means

% Chose to compute the effect size over pairs (difference against 0)

% plots stored means and SE since they reflect the general anova that was done, 
% so for table values here you should still make CIs

save_file = ['chance_prop_' general_save_file]

rois_col = [];
Ns = [];
Zs = [];
rs = [];
ps = [];
medians = [];
chance_medians = [];
chance_props = [];

for run = 1:length(map_types)
    ROIs = rois.(map_types{run});
    if string(map_types{run}) == "num"
        DT_run = 'NumerosityAll';
    else
        DT_run = 'TimingAll';
    end

    for roi = 1:length(ROIs)
        rois_col = [rois_col; {ROIs{roi}}];

        data_this_roi = stat.data.(DT_run).prop(string(stat.data.(DT_run).map) == ROIs{roi});

        Ns = [Ns; {num2str(sum(~isnan(data_this_roi)))}];

        median = nanmedian(data_this_roi);
        % make CIs
        bootstr_CI = bootci(1000, {@nanmedian, data_this_roi},'alpha',0.05);
        CI_prop = [num2str(median,'%.2f'), ' [', num2str(bootstr_CI(1),'%.2f'), ', ' num2str(bootstr_CI(2),'%.2f'), ']'];
        medians = [medians; {CI_prop}];


        chance_data_this_roi = stat.data.(DT_run).chance_proportion(string(stat.data.(DT_run).map) == ROIs{roi});
        chance_median = nanmedian(chance_data_this_roi);
        chance_bootstr_CI = bootci(1000, {@nanmedian, chance_data_this_roi},'alpha',0.05);
        chance_CI_prop = [num2str(chance_median,'%.2f'), ' [', num2str(chance_bootstr_CI(1),'%.2f'), ', ' num2str(chance_bootstr_CI(2),'%.2f'), ']'];
        chance_medians = [chance_medians; {chance_CI_prop}];


        Zval = num2str(stat.different_from_chance.(DT_run).proportion.proportion.(ROIs{roi}).stats.zval,'%.2f');

        Zs = [Zs; {Zval}];

        N = sum(~isnan(stat.data.(DT_run).prop(stat.data.(DT_run).map==ROIs{roi})));
        eff_size = num2str(stat.different_from_chance.(DT_run).proportion.proportion.(ROIs{roi}).stats.zval/sqrt(N),'%.2f');

        rs = [rs; {eff_size}];
        p = num2str(stat.different_from_chance.(DT_run).proportion.proportion.(ROIs{roi}).adj_p,'%.3f');

        ps = [ps; {p}];

    end
end

tbl = table();
tbl.rois = categorical(rois_col);
tbl.N = categorical(Ns);
tbl.median_CI = categorical(medians);
tbl.chance_median_CI = categorical(chance_medians);
tbl.Z = categorical(Zs);
tbl.r = categorical(rs);
tbl.p = categorical(ps);

writetable(tbl, save_file)

%% post hocs

save_file = ['posthoc_prop_' general_save_file]

rois_1_col = [];
rois_2_col = [];
diffs = [];
ps = [];
means_roi_1 = [];
means_roi_2 = [];

for run = 1:length(map_types)
    ROIs = rois.(map_types{run});
    if string(map_types{run}) == "num"
        DT_run = 'NumerosityAll';
    else
        DT_run = 'TimingAll';
    end

    for roi_1 = 1:length(ROIs)-1
        for roi_2 = roi_1 + 1: length(ROIs)

            rois_1_col = [rois_1_col; {ROIs{roi_1}}];
            rois_2_col = [rois_2_col; {ROIs{roi_2}}];

            mean_1 = plot_values.(DT_run).proportion.point(roi_1);
            mean_2 = plot_values.(DT_run).proportion.point(roi_2);

            mean_SE_1 = [mean_1-plot_values.(DT_run).proportion.errorbars(1,roi_1) mean_1+plot_values.(DT_run).proportion.errorbars(1,roi_1)];
            mean_SE_2 = [mean_2-plot_values.(DT_run).proportion.errorbars(1,roi_2) mean_2+plot_values.(DT_run).proportion.errorbars(1,roi_2)];

            CI_1 = [num2str(mean_1,'%.2f'), ' [', num2str(mean_SE_1(1),'%.2f'), ', ' num2str(mean_SE_1(2),'%.2f'), ']'];
            CI_2 = [num2str(mean_2,'%.2f'), ' [', num2str(mean_SE_2(1),'%.2f'), ', ' num2str(mean_SE_2(2),'%.2f'), ']'];

            means_roi_1 = [means_roi_1; {CI_1}];
            means_roi_2 = [means_roi_2; {CI_2}];

            roi_id_1 = find(contains(stat.posthoc.(DT_run).proportion.gnames,ROIs{roi_1}));
            roi_id_2 = find(contains(stat.posthoc.(DT_run).proportion.gnames,ROIs{roi_2}));

            if roi_id_1 > roi_id_2
                d_diff = -stat.posthoc.(DT_run).proportion.multcomp(stat.posthoc.(DT_run).proportion.multcomp(:,1) == roi_id_2 & stat.posthoc.(DT_run).proportion.multcomp(:,2) == roi_id_1,4);
                ll_diff = -stat.posthoc.(DT_run).proportion.multcomp(stat.posthoc.(DT_run).proportion.multcomp(:,1) == roi_id_2 & stat.posthoc.(DT_run).proportion.multcomp(:,2) == roi_id_1,3);
                ul_diff = -stat.posthoc.(DT_run).proportion.multcomp(stat.posthoc.(DT_run).proportion.multcomp(:,1) == roi_id_2 & stat.posthoc.(DT_run).proportion.multcomp(:,2) == roi_id_1,5);

                p = num2str(stat.posthoc.(DT_run).proportion.multcomp(stat.posthoc.(DT_run).proportion.multcomp(:,1) == roi_id_2 & stat.posthoc.(DT_run).proportion.multcomp(:,2) == roi_id_1,6),'%.3f');

            else

                d_diff = stat.posthoc.(DT_run).proportion.multcomp(stat.posthoc.(DT_run).proportion.multcomp(:,1) == roi_id_1 & stat.posthoc.(DT_run).proportion.multcomp(:,2) == roi_id_2,4);
                ll_diff = stat.posthoc.(DT_run).proportion.multcomp(stat.posthoc.(DT_run).proportion.multcomp(:,1) == roi_id_1 & stat.posthoc.(DT_run).proportion.multcomp(:,2) == roi_id_2,3);
                ul_diff = stat.posthoc.(DT_run).proportion.multcomp(stat.posthoc.(DT_run).proportion.multcomp(:,1) == roi_id_1 & stat.posthoc.(DT_run).proportion.multcomp(:,2) == roi_id_2,5);

                p = num2str(stat.posthoc.(DT_run).proportion.multcomp(stat.posthoc.(DT_run).proportion.multcomp(:,1) == roi_id_1 & stat.posthoc.(DT_run).proportion.multcomp(:,2) == roi_id_2,6),'%.3f');

            end
            diff =  [num2str(d_diff,'%.2f'), ' [', num2str(ll_diff,'%.2f'), ', ' num2str(ul_diff,'%.2f'), ']'];
            diffs = [diffs; {diff}];



            ps = [ps; {p}];
        end

    end
end

tbl = table();
tbl.roi_1 = categorical(rois_1_col);
tbl.roi_2 = categorical(rois_2_col);
tbl.means_1 = categorical(means_roi_1);
tbl.means_2 = categorical(means_roi_2);
tbl.difference_in_mutlcompare = categorical(diffs);
tbl.p = categorical(ps);

writetable(tbl, save_file)

%% repeated measures
clear all

whichCombi = '4'
file_general = ['NumerosityAll_TimingAll_minVEtime=0.2_minVEnum=0.3_whichCombi=' whichCombi '_topo_measure=mean'];
general_save_file = ['data_tbl_' file_general '.csv'];

map_types = {'num','time'};

% rois
rois.time = ["TLO","TTOP","TTOA","TPO","TLS","TPCI","TPCM","TPCS","TFI","TFS"];
rois.num = ["NLO","NTO","NPO","NPCI","NPCM","NPCS","NFI","NFS"];

% load the rm files
load('compare_rm_proportion_minVEtime=0.2_minVEnum=0.3_whichCombi=4_stat.mat')
load('compare_rm_proportion_minVEtime=0.2_minVEnum=0.3_whichCombi=4_plot_values.mat')

save_file = ['rm_prop_' general_save_file]

rois_col = [];
ps = [];
medians_roi_1 = [];
medians_roi_2 = [];
Zs = [];
Ns = [];
rs = [];
ps = [];

for run = 1:length(map_types)
    ROIs = rois.(map_types{run});
    if string(map_types{run}) == "num"
        DT_run = 'NumerosityAll';
    else
        DT_run = 'TimingAll';
    end

    for roi = 1:length(ROIs)

        rois_col = [rois_col; {ROIs{roi}}];

        median_1 = plot_values.(map_types{run}).proportion.point.diff(roi);
        median_2 = plot_values.(map_types{run}).proportion.point.same(roi);

        median_CI_1 = [median_1-plot_values.(map_types{run}).proportion.errorbars.diff(1,roi) median_1+plot_values.(map_types{run}).proportion.errorbars.diff(2,roi)];
        median_CI_2 = [median_2-plot_values.(map_types{run}).proportion.errorbars.same(1,roi) median_2+plot_values.(map_types{run}).proportion.errorbars.same(2,roi)];

        CI_1 = [num2str(median_1,'%.2f'), ' [', num2str(median_CI_1(1),'%.2f'), ', ' num2str(median_CI_1(2),'%.2f'), ']'];
        CI_2 = [num2str(median_2,'%.2f'), ' [', num2str(median_CI_2(1),'%.2f'), ', ' num2str(median_CI_2(2),'%.2f'), ']'];

        medians_roi_1 = [medians_roi_1; {CI_1}];
        medians_roi_2 = [medians_roi_2; {CI_2}];


        Zval = num2str(stat.(map_types{run}).proportion.(ROIs{roi}).stats.zval,'%.2f');
        Zs = [Zs; {Zval}];

        N = sum(~isnan(plot_values.(map_types{run}).proportion.scatter.prop.diff(plot_values.(map_types{run}).proportion.scatter.maps==ROIs{roi})) & ~isnan(plot_values.(map_types{run}).proportion.scatter.prop.same(plot_values.(map_types{run}).proportion.scatter.maps==ROIs{roi})));
        eff_size = num2str(stat.(map_types{run}).proportion.(ROIs{roi}).stats.zval/sqrt(N),'%.2f');

        Ns = [Ns; {num2str(N)}]

        rs = [rs; {eff_size}];
        p = num2str(stat.(map_types{run}).proportion.(ROIs{roi}).adj_p,'%.3f');

        ps = [ps; {p}];


    end
end

tbl = table();
tbl.roi = categorical(rois_col);
tbl.N = categorical(Ns);
tbl.medians_num_time = categorical(medians_roi_1);
tbl.medians_rm = categorical(medians_roi_2);
tbl.Z = categorical(Zs);
tbl.r = categorical(rs);
tbl.p = categorical(ps);


writetable(tbl, save_file)



%% correlation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
whichCombi = '4';
file_general = ['NumerosityAll_TimingAll_minVEtime=0.2_minVEnum=0.3_whichCombi=' whichCombi '_topo_measure=mean'];
stat_file = ['stat_' file_general '.mat'];
general_save_file = ['data_tbl_' file_general '.csv'];

load(stat_file)
% load plot values (including indicated means and medians)
load(['minVEtime=0.2_minVEnum=0.3_whichCombi=',whichCombi,'_correlation_NumerosityAll_TimingAll_plot_values.mat'])


map_types = {'num','time'};

%% rois
rois.time = ["TLO","TTOP","TTOA","TPO","TLS","TPCI","TPCM","TPCS","TFI","TFS"];
rois.num = ["NLO","NTO","NPO","NPCI","NPCM","NPCS","NFI","NFS"];


%% different from chance
save_file = ['chance_corr_' general_save_file]

rois_col = [{'Numerosity vs Duration'}];
Ns = [{''}];
Zs = [{''}];
rs = [{''}];
ps = [{''}];
medians = [{''}];

parameters = {'x0x0','x0y0'};

for par = 1:length(parameters)
    for run = 1:length(map_types)
        ROIs = rois.(map_types{run});
        if string(map_types{run}) == "num"
            DT_run = 'NumerosityAll';
        else
            DT_run = 'TimingAll';
        end


        for roi = 1:length(ROIs)
            data_this_roi = stat.data.(DT_run).correlation.x0.(parameters{par}(end-1:end))(string(stat.data.(DT_run).map) == ROIs{roi});
            N = sum(~isnan(data_this_roi));

            rois_col = [rois_col; {ROIs{roi}}];
            Ns = [Ns; {num2str(N)}];

            if N <= 1
                medians = [medians; {num2str(data_this_roi(~isnan(data_this_roi)),'%.2f')}];
                Zs = [Zs; {''}];
                rs = [rs; {''}];
                ps = [ps; {''}];
                continue
            end

            if strcmp(whichCombi, '4') == 1
                median = plot_values.(DT_run).correlation.(parameters{par}).point(roi);
                median_CI = [median-plot_values.(DT_run).correlation.(parameters{par}).errorbars(1,roi) median+plot_values.(DT_run).correlation.(parameters{par}).errorbars(2,roi)];
                CI = [num2str(median,'%.2f'), ' [', num2str(median_CI(1),'%.2f'), ', ' num2str(median_CI(2),'%.2f'), ']'];
                medians = [medians; {CI}];
            elseif strcmp(whichCombi, '2') == 1
                median = nanmedian(data_this_roi);
                bootstr_CI = bootci(1000, {@nanmedian, data_this_roi},'alpha',0.05);
                CI_corr = [num2str(median,'%.2f'), ' [', num2str(bootstr_CI(1),'%.2f'), ', ' num2str(bootstr_CI(2),'%.2f'), ']'];
                medians = [medians; {CI_corr}];
            end

            Zval = num2str(stat.different_from_chance.(DT_run).correlation.(parameters{par}).(ROIs{roi}).stats.zval,'%.2f');
            Zs = [Zs; {Zval}];

            eff_size = num2str(stat.different_from_chance.(DT_run).correlation.(parameters{par}).(ROIs{roi}).stats.zval/sqrt(N),'%.2f');
            rs = [rs; {eff_size}];
            p = num2str(stat.different_from_chance.(DT_run).correlation.(parameters{par}).(ROIs{roi}).adj_p,'%.3f');

            ps = [ps; {p}];

        end
    end

    if par == 1
        rois_col = [rois_col; {'Numerosity vs Period'}];
        Ns = [Ns;  {''}];
        medians = [medians; {''}];
        Zs = [Zs; {''}];
        rs = [rs; {''}];
        ps = [ps; {''}];
    end
end


tbl = table();
tbl.rois = categorical(rois_col);
tbl.N = categorical(Ns);
tbl.correlation_median = categorical(medians);
tbl.Z = categorical(Zs);
tbl.r = categorical(rs);
tbl.p = categorical(ps);

writetable(tbl, save_file)


%% post hocs whichCombi = 4
clear all
whichCombi = '4';
file_general = ['NumerosityAll_TimingAll_minVEtime=0.2_minVEnum=0.3_whichCombi=' whichCombi '_topo_measure=mean'];
stat_file = ['stat_' file_general '.mat'];
general_save_file = ['data_tbl_' file_general '.csv'];

load(stat_file)
load('minVEtime=0.2_minVEnum=0.3_whichCombi=4_correlation_NumerosityAll_TimingAll_plot_values.mat')


map_types = {'num','time'};

% rois
rois.time = ["TLO","TTOP","TTOA","TPO","TLS","TPCI","TPCM","TPCS","TFI","TFS"];
rois.num = ["NLO","NTO","NPO","NPCI","NPCM","NPCS","NFI","NFS"];


save_file = ['posthoc_corr_' general_save_file]

rois_col_1 = [];
medians_1 = [];
rois_col_2 = [];
medians_2 = [];
Qs = [];


parameters = {'x0x0','x0y0'};

for run = 1:length(map_types)
    for par = 1:length(parameters)
        if par == 1
            rois_col_1 = [rois_col_1;{'Numerosity vs Duration'}];
            medians_1 = [medians_1;{''}];
            rois_col_2 = [rois_col_2;{''}];
            medians_2 = [medians_2;{''}];
            Qs = [Qs;{''}];
        end


        ROIs = rois.(map_types{run});
        if string(map_types{run}) == "num"
            DT_run = 'NumerosityAll';
        else
            DT_run = 'TimingAll';
        end

        Dunn_tbl = stat.posthoc.(DT_run).correlation.(parameters{par}).tbl;
        Dunn_rois = stat.posthoc.(DT_run).correlation.(parameters{par}).roi_names;

        size_Dunn = size(Dunn_tbl)

        for roi_rows = 1:size_Dunn(1)

            roi_ids = split(Dunn_tbl{roi_rows,1},'-');
            roi_id_1 = str2num(roi_ids{1});
            roi_id_2 = str2num(roi_ids{2});


            rois_col_1 = [rois_col_1; {Dunn_rois{roi_id_1}}];
            rois_col_2 = [rois_col_2; {Dunn_rois{roi_id_2}}];
            [~,~,roi_id_1_plot_values] = intersect(Dunn_rois{roi_id_1},rois.(map_types{run}));
            [~,~,roi_id_2_plot_values] = intersect(Dunn_rois{roi_id_2},rois.(map_types{run}));

            median_1 = plot_values.(DT_run).correlation.(parameters{par}).point(roi_id_1_plot_values);
            median_CI_1 = [median_1-plot_values.(DT_run).correlation.(parameters{par}).errorbars(1,roi_id_1_plot_values) median_1+plot_values.(DT_run).correlation.(parameters{par}).errorbars(2,roi_id_1_plot_values)];
            CI_1 = [num2str(median_1,'%.2f'), ' [', num2str(median_CI_1(1),'%.2f'), ', ' num2str(median_CI_1(2),'%.2f'), ']'];
            medians_1 = [medians_1; {CI_1}];

            median_2 = plot_values.(DT_run).correlation.(parameters{par}).point(roi_id_2_plot_values);
            median_CI_2 = [median_2-plot_values.(DT_run).correlation.(parameters{par}).errorbars(1,roi_id_2_plot_values) median_2+plot_values.(DT_run).correlation.(parameters{par}).errorbars(2,roi_id_2_plot_values)];
            CI_2 = [num2str(median_2,'%.2f'), ' [', num2str(median_CI_2(1),'%.2f'), ', ' num2str(median_CI_2(2),'%.2f'), ']'];
            medians_2 = [medians_2; {CI_2}];

            Qval = num2str(Dunn_tbl{roi_rows,2},'%.2f');
            Qs = [Qs; {Qval}];

        end


        if par == 1
            rois_col_1 = [rois_col_1;{'Numerosity vs Period'}];
            medians_1 = [medians_1;{''}];
            rois_col_2 = [rois_col_2;{''}];
            medians_2 = [medians_2;{''}];
            Qs = [Qs;{''}];
        end
    end
end


tbl = table();
tbl.rois_1 = categorical(rois_col_1);
tbl.correlation_median_1 = categorical(medians_1);
tbl.rois_2 = categorical(rois_col_2);
tbl.correlation_median_2 = categorical(medians_2);
tbl.Q = categorical(Qs);

writetable(tbl, save_file)

%% post hoc correlation which combi = 2

clear all
whichCombi = '2';
file_general = ['NumerosityAll_TimingAll_minVEtime=0.2_minVEnum=0.3_whichCombi=' whichCombi '_topo_measure=mean'];
stat_file = ['stat_' file_general '.mat'];
general_save_file = ['data_tbl_' file_general '.csv'];

load(stat_file)
load('minVEtime=0.2_minVEnum=0.3_whichCombi=2_correlation_NumerosityAll_TimingAll_plot_values.mat')


map_types = {'num','time'};

% rois
rois.time = ["TLO","TTOP","TTOA","TPO","TLS","TPCI","TPCM","TPCS","TFI","TFS"];
rois.num = ["NLO","NTO","NPO","NPCI","NPCM","NPCS","NFI","NFS"];

save_file = ['posthoc_corr_' general_save_file];

rois_1_col = [];
rois_2_col = [];
diffs = [];
ps = [];
means_roi_1 = [];
means_roi_2 = [];

parameters = {'x0x0','x0y0'};

for run = 1:length(map_types)
    ROIs = rois.(map_types{run});
    if string(map_types{run}) == "num"
        DT_run = 'NumerosityAll';
    else
        DT_run = 'TimingAll';
    end

    for par = 1:length(parameters)

        if par == 1
            rois_1_col = [rois_1_col;{'Numerosity vs Duration'}];
            means_roi_1 = [means_roi_1;{''}];
            rois_2_col = [rois_2_col;{''}];
            means_roi_2 = [means_roi_2;{''}];
            diffs = [diffs;{''}];
            ps = [ps;{''}];
        end

        for roi_1 = 1:length(ROIs)-1
            for roi_2 = roi_1 + 1: length(ROIs)

                rois_1_col = [rois_1_col; {ROIs{roi_1}}];
                rois_2_col = [rois_2_col; {ROIs{roi_2}}];

                mean_1 = plot_values.(DT_run).correlation.(parameters{par}).point(roi_1);
                mean_2 = plot_values.(DT_run).correlation.(parameters{par}).point(roi_2);

                mean_SE_1 = [mean_1-plot_values.(DT_run).correlation.(parameters{par}).errorbars(1,roi_1) mean_1+plot_values.(DT_run).correlation.(parameters{par}).errorbars(1,roi_1)];
                mean_SE_2 = [mean_2-plot_values.(DT_run).correlation.(parameters{par}).errorbars(1,roi_2) mean_2+plot_values.(DT_run).correlation.(parameters{par}).errorbars(1,roi_2)];

                CI_1 = [num2str(mean_1,'%.2f'), ' [', num2str(mean_SE_1(1),'%.2f'), ', ' num2str(mean_SE_1(2),'%.2f'), ']'];
                CI_2 = [num2str(mean_2,'%.2f'), ' [', num2str(mean_SE_2(1),'%.2f'), ', ' num2str(mean_SE_2(2),'%.2f'), ']'];

                means_roi_1 = [means_roi_1; {CI_1}];
                means_roi_2 = [means_roi_2; {CI_2}];

                roi_id_1 = find(contains(stat.posthoc.(DT_run).correlation.(parameters{par}).gnames,ROIs{roi_1}));
                roi_id_2 = find(contains(stat.posthoc.(DT_run).correlation.(parameters{par}).gnames,ROIs{roi_2}));

                if isempty(roi_id_1) || isempty(roi_id_2)
                    diffs = [diffs; {''}];



                    ps = [ps; {''}];
                    continue
                end

                if roi_id_1 > roi_id_2
                    d_diff = -stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,1) == roi_id_2 & stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,2) == roi_id_1,4);
                    ll_diff = -stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,1) == roi_id_2 & stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,2) == roi_id_1,3);
                    ul_diff = -stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,1) == roi_id_2 & stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,2) == roi_id_1,5);

                    p = num2str(stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,1) == roi_id_2 & stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,2) == roi_id_1,6),'%.3f');

                else

                    d_diff = stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,1) == roi_id_1 & stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,2) == roi_id_2,4);
                    ll_diff = stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,1) == roi_id_1 & stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,2) == roi_id_2,3);
                    ul_diff = stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,1) == roi_id_1 & stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,2) == roi_id_2,5);

                    p = num2str(stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,1) == roi_id_1 & stat.posthoc.(DT_run).correlation.(parameters{par}).multcomp(:,2) == roi_id_2,6),'%.3f');

                end
                diff =  [num2str(d_diff,'%.2f'), ' [', num2str(ll_diff,'%.2f'), ', ' num2str(ul_diff,'%.2f'), ']'];
                diffs = [diffs; {diff}];



                ps = [ps; {p}];
            end

        end
        if par == 1
            rois_1_col = [rois_1_col;{'Numerosity vs Period'}];
            means_roi_1 = [means_roi_1;{''}];
            rois_2_col = [rois_2_col;{''}];
            means_roi_2 = [means_roi_2;{''}];
            diffs = [diffs;{''}];
            ps = [ps;{''}];
        end
    end
end

tbl = table();
tbl.roi_1 = categorical(rois_1_col);
tbl.roi_2 = categorical(rois_2_col);
tbl.means_1 = categorical(means_roi_1);
tbl.means_2 = categorical(means_roi_2);
tbl.difference_in_mutlcompare = categorical(diffs);
tbl.p = categorical(ps);

writetable(tbl, save_file)


%% repeated measures chance
clear all

whichCombi = '5'
file_general = ['NumerosityAll_TimingAll_minVEtime=0.2_minVEnum=0.3_whichCombi=' whichCombi '_topo_measure=mean'];
general_save_file = ['data_tbl_' file_general '.csv'];

map_types = {'num','time'};

% rois
rois.time = ["TLO","TTOP","TTOA","TPO","TLS","TPCI","TPCM","TPCS","TFI","TFS"];
rois.num = ["NLO","NTO","NPO","NPCI","NPCM","NPCS","NFI","NFS"];


load(['compare_rm_correlation_minVEtime=0.2_minVEnum=0.3_whichCombi=',whichCombi,'_stat.mat'])
load(['compare_rm_correlation_minVEtime=0.2_minVEnum=0.3_whichCombi=',whichCombi,'_plot_values.mat'])

save_file = ['rm_corr_chance_' general_save_file]

rois_col = [{'Numerosity vs Duration'}];
Ns = [{''}];
Zs = [{''}];
rs = [{''}];
ps = [{''}];
medians = [{''}];

diff_parameters = {'x0x0','x0y0'};
same_parameters = {'x0x0','y0y0'}

for par = 1:length(diff_parameters)
    for run = 1:length(map_types)
        if run == 1 && par == 2
            continue
        end


        ROIs = rois.(map_types{run});
        if string(map_types{run}) == "num"
            DT_run = 'NumerosityAll';
        else
            DT_run = 'TimingAll';
        end

        roi_id = 0;
        for roi = 1:length(ROIs)


            data_this_roi_same = plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{par}).scatter.corr.same(string(plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{par}).scatter.maps) == ROIs{roi});
            data_this_roi_diff = plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{par}).scatter.corr.diff(string(plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{par}).scatter.maps) == ROIs{roi});

            data_this_roi_same(isnan(data_this_roi_diff)) = NaN;
            data_this_roi_diff(isnan(data_this_roi_same)) = NaN;

            N = sum(~isnan(data_this_roi_same));



            if N <= 1
                continue
            end

            roi_id = roi_id + 1;
            rois_col = [rois_col; {ROIs{roi}}];
            Ns = [Ns; {num2str(N)}];

            median = plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{par}).point.same(roi_id);
            median_CI = [median-plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{par}).errorbars.same(1,roi_id) median+plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{par}).errorbars.same(2,roi_id)];
            CI = [num2str(median,'%.2f'), ' [', num2str(median_CI(1),'%.2f'), ', ' num2str(median_CI(2),'%.2f'), ']'];
            medians = [medians; {CI}];


            Zval = num2str(stat.(map_types{run}).chance_distribution.correlation.(same_parameters{par}).(ROIs{roi}).stats.zval,'%.2f');
            Zs = [Zs; {Zval}];

            eff_size = num2str(stat.(map_types{run}).chance_distribution.correlation.(same_parameters{par}).(ROIs{roi}).stats.zval/sqrt(N),'%.2f');
            rs = [rs; {eff_size}];
            p = num2str(stat.(map_types{run}).chance_distribution.correlation.(same_parameters{par}).(ROIs{roi}).adj_p,'%.3f');

            ps = [ps; {p}];

        end
    end

    if par == 1
        rois_col = [rois_col; {'Numerosity vs Period'}];
        Ns = [Ns;  {''}];
        medians = [medians; {''}];
        Zs = [Zs; {''}];
        rs = [rs; {''}];
        ps = [ps; {''}];
    end
end


tbl = table();
tbl.rois = categorical(rois_col);
tbl.N = categorical(Ns);
tbl.correlation_median = categorical(medians);
tbl.Z = categorical(Zs);
tbl.r = categorical(rs);
tbl.p = categorical(ps);

writetable(tbl, save_file)

%% repeated measures pairwise

clear all

whichCombi = '5'
file_general = ['NumerosityAll_TimingAll_minVEtime=0.2_minVEnum=0.3_whichCombi=' whichCombi '_topo_measure=mean'];
general_save_file = ['data_tbl_' file_general '.csv'];

map_types = {'num','time'};

% rois
rois.time = ["TLO","TTOP","TTOA","TPO","TLS","TPCI","TPCM","TPCS","TFI","TFS"];
rois.num = ["NLO","NTO","NPO","NPCI","NPCM","NPCS","NFI","NFS"];


load(['compare_rm_correlation_minVEtime=0.2_minVEnum=0.3_whichCombi=',whichCombi,'_stat.mat'])
load(['compare_rm_correlation_minVEtime=0.2_minVEnum=0.3_whichCombi=',whichCombi,'_plot_values.mat'])

save_file = ['rm_corr_' general_save_file]

rois_col = [];
ps = [];
medians_roi_same = [];
medians_roi_diff = [];
Zs = [];
Ns = [];
rs = [];
ps = [];

diff_parameters = {'x0x0','x0y0'};
same_parameters = {'x0x0','y0y0'};


for par = 1:length(diff_parameters)
    for run = 1:length(map_types)
        if par == 1 && run == 1
            rois_col = [rois_col; {'Numerosity-Duration vs Numerosity-Numerosity'}];
            Ns = [Ns;  {''}];
            medians_roi_diff = [medians_roi_diff; {''}];
            medians_roi_same = [medians_roi_same; {''}];
            Zs = [Zs; {''}];
            rs = [rs; {''}];
            ps = [ps; {''}];
        end



        if run == 2 && par == 2
            spar = 2;
        else
            spar = 1;
        end


        ROIs = rois.(map_types{run});
        if string(map_types{run}) == "num"
            DT_run = 'NumerosityAll';
        else
            DT_run = 'TimingAll';
        end

        roi_id = 0;
        for roi = 1:length(ROIs)


            data_this_roi_same = plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{spar}).scatter.corr.same(string(plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{spar}).scatter.maps) == ROIs{roi});
            data_this_roi_diff = plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{spar}).scatter.corr.diff(string(plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{spar}).scatter.maps) == ROIs{roi});

            data_this_roi_same(isnan(data_this_roi_diff)) = NaN;
            data_this_roi_diff(isnan(data_this_roi_same)) = NaN;

            N = sum(~isnan(data_this_roi_same));


            if N <= 1
                continue
            end


            roi_id = roi_id + 1;
            rois_col = [rois_col; {ROIs{roi}}];
            Ns = [Ns; {num2str(N)}];

            median_same = plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{spar}).point.same(roi_id);
            median_CI_same = [median_same-plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{spar}).errorbars.same(1,roi_id) median_same+plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{spar}).errorbars.same(2,roi_id)];
            CI_same = [num2str(median_same,'%.2f'), ' [', num2str(median_CI_same(1),'%.2f'), ', ' num2str(median_CI_same(2),'%.2f'), ']'];
            medians_roi_same = [medians_roi_same; {CI_same}];

            median_diff = plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{spar}).point.diff(roi_id);
            median_CI_diff = [median_diff-plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{spar}).errorbars.diff(1,roi_id) median_diff+plot_values.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{spar}).errorbars.diff(2,roi_id)];
            CI_diff = [num2str(median_diff,'%.2f'), ' [', num2str(median_CI_diff(1),'%.2f'), ', ' num2str(median_CI_diff(2),'%.2f'), ']'];
            medians_roi_diff = [medians_roi_diff; {CI_diff}];

            %%
            Zval = num2str(stat.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{spar}).(ROIs{roi}).stats.zval,'%.2f');
            Zs = [Zs; {Zval}];

            eff_size = num2str(stat.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{spar}).(ROIs{roi}).stats.zval/sqrt(N),'%.2f');
            rs = [rs; {eff_size}];


            p = num2str(stat.(map_types{run}).correlation.(diff_parameters{par}).(same_parameters{spar}).(ROIs{roi}).adj_p,'%.3f');
            ps = [ps; {p}];



        end

        if par == 1 && run == 1
            rois_col = [rois_col; {'Numerosity-Duration vs Duration-Duration'}];
            Ns = [Ns;  {''}];
            medians_roi_diff = [medians_roi_diff; {''}];
            medians_roi_same = [medians_roi_same; {''}];
            Zs = [Zs; {''}];
            rs = [rs; {''}];
            ps = [ps; {''}];
        end

        if par == 1 && run == 2
            rois_col = [rois_col; {'Numerosity-Period vs Numerosity-Numerosity'}];
            Ns = [Ns;  {''}];
            medians_roi_diff = [medians_roi_diff; {''}];
            medians_roi_same = [medians_roi_same; {''}];
            Zs = [Zs; {''}];
            rs = [rs; {''}];
            ps = [ps; {''}];
        end

        if par == 2 && run == 1
            rois_col = [rois_col; {'Numerosity-Period vs Period-Period'}];
            Ns = [Ns;  {''}];
            medians_roi_diff = [medians_roi_diff; {''}];
            medians_roi_same = [medians_roi_same; {''}];
            Zs = [Zs; {''}];
            rs = [rs; {''}];
            ps = [ps; {''}];
        end
    end
end


tbl = table();
tbl.roi = categorical(rois_col);
tbl.N = categorical(Ns);
tbl.medians_num_time = categorical(medians_roi_diff);
tbl.medians_rm = categorical(medians_roi_same);
tbl.Z = categorical(Zs);
tbl.r = categorical(rs);
tbl.p = categorical(ps);


writetable(tbl, save_file)


%% TOPO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
whichCombi = '4';
DT_runs = {'NumerosityAll','TimingAll'};

if contains(DT_runs{1},"Numerosity") &&  contains(DT_runs{2},"Numerosity")
     csv_addition = '_rm_numerosity';
     parameters = {'x0x0'};

elseif contains(DT_runs{1},"Timing") &&  contains(DT_runs{2},"Timing")
     csv_addition = '_rm_timing';
     parameters = {'x0x0','y0y0'};
else
    csv_addition = '';
    parameters = {'x0x0','x0y0'};
end
csv_addition = '_leftRight';


file_general = [DT_runs{1}, '_', DT_runs{2}, '_minVEtime=0.2_minVEnum=0.3_whichCombi=', whichCombi, '_topo_measure=mean'];
% % stat_file = ['stat_' file_general '.mat'];
stat_file = ['left-right-topo_stat_' file_general '.mat'];
general_save_file = ['data_tbl_' file_general '.csv'];

load(stat_file)
% % load(['minVEtime=0.2_minVEnum=0.3_whichCombi=',whichCombi,'_topo_', DT_runs{1}, '_', DT_runs{2},'_plot_values.mat'])
load(['left-right-topo_minVEtime=0.2_minVEnum=0.3_whichCombi=',whichCombi,'_topo_', DT_runs{1}, '_', DT_runs{2},'_plot_values.mat'])

%% rois
rois.time = ["TLO","TTOP","TTOA","TPO","TLS","TPCI","TPCM","TPCS","TFI","TFS"];
rois.num = ["NLO","NTO","NPO","NPCI","NPCM","NPCS","NFI","NFS"];


%% different from chance
% % save_file = ['BF_topo_' general_save_file];
save_file = ['BF_left_right_topo_' general_save_file];

rois_col = [];
Ns = [];
BFs = [];
H0s = [];
HAs = [];
means = [];

for par = 1:length(parameters)
    rois_col = [rois_col;{parameters{par}}];
    Ns = [Ns;{''}];
    BFs = [BFs;{''}];
    H0s = [H0s;{''}];
    HAs = [HAs;{''}];
    means = [means;{''}];
    for run = 1:length(DT_runs)
        DT_run = DT_runs{run};

        if contains(DT_runs{1},"Numerosity") &&  contains(DT_runs{2},"Numerosity")
            if run == 2
                continue
            end
            DT_run = 'NumerosityHalvesMean';

        elseif contains(DT_runs{1},"Timing") &&  contains(DT_runs{2},"Timing")
            if run == 2
                continue
            end
            DT_run = 'TimingHalvesMean';
        end

        if contains(DT_run,"Numerosity") 
            ROIs = rois.num;
        else
            ROIs = rois.time;
        end

        for roi = 1:length(ROIs)

            if (contains(DT_runs{1},"Numerosity") &&  contains(DT_runs{2},"Numerosity")) || (contains(DT_runs{1},"Timing") &&  contains(DT_runs{2},"Timing"))
                % only used for N and if 1 datapoint mean
                data_half = stat.data.(DT_runs{1}).topo.(parameters{par}(1:2)).(parameters{par}(3:4))(string(stat.data.(DT_runs{1}).map) == ROIs{roi});
                data_other_half = stat.data.(DT_runs{2}).topo.(parameters{par}(1:2)).(parameters{par}(3:4))(string(stat.data.(DT_runs{2}).map) == ROIs{roi});
                data_half(isnan(data_half)&~isnan(data_other_half))=data_other_half(isnan(data_half)&~isnan(data_other_half));
                data_other_half(~isnan(data_half)&isnan(data_other_half)) = data_half(~isnan(data_half)&isnan(data_other_half));

                data_this_roi = rad2deg(circ_mean([deg2rad(data_half);deg2rad(data_other_half)]',[],2));

            else
                data_this_roi = stat.data.(DT_run).topo.(parameters{par}(1:2)).(parameters{par}(3:4))(string(stat.data.(DT_run).map) == ROIs{roi});
            end

            N = sum(~isnan(data_this_roi));

            rois_col = [rois_col; {ROIs{roi}}];
            Ns = [Ns; {num2str(N)}];

            if N <= 1
                mean = data_this_roi(~isnan(data_this_roi));
                if mean > 180
                    mean = (360 - mean)*-1;
                end

                means = [means; {num2str(mean,'%.2f')}];
                BFs = [BFs;{''}];
                H0s = [H0s;{''}];
                HAs = [HAs;{''}];
                continue
            end

                mean = plot_values.(DT_run).topo.(parameters{par}).(ROIs{roi}).mean;
                if mean > 180
                    mean = (360 - mean)*-1;
                end

                SE_1 = plot_values.(DT_run).topo.(parameters{par}).(ROIs{roi}).errorbars(1);
                if SE_1 > 180
                    SE_1 = (360 - SE_1)*-1;
                end

                SE_2 = plot_values.(DT_run).topo.(parameters{par}).(ROIs{roi}).errorbars(2);
                if SE_2 > 180
                    SE_2 = (360 - SE_2)*-1;
                end

                CI = [num2str(mean,'%.2f'), ' [', num2str(SE_1,'%.2f'), ', ' num2str(SE_2,'%.2f'), ']'];
                means = [means; {CI}];

                % Bayes factor info
                csv_file = ['BF_wc', whichCombi, '_', parameters{par}];

                BF_info = readcell([csv_file csv_addition '.csv']);
                
                roi_id = find(BF_info(:,2) == ROIs(roi));
                pH0 = num2str(BF_info{roi_id,3},'%.2f');
                pHA = num2str(BF_info{roi_id,4},'%.2f');
                BF = num2str(BF_info{roi_id,5},'%.3f');

                H0s = [H0s; {pH0}];
                HAs = [HAs; {pHA}];
                BFs = [BFs; {BF}];
              
        end
    end
end


tbl = table();
tbl.rois = categorical(rois_col);
tbl.N = categorical(Ns);
tbl.angle_means = categorical(means);
tbl.pH0 = categorical(H0s);
tbl.pHa = categorical(HAs);
tbl.BFs = categorical(BFs);

writetable(tbl, save_file)
