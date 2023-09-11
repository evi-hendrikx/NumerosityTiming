function NumTim_compare_correlations(fig_path, whichCombi)
%% Compares correlations for Num-Time vs repeated measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% fig_path: path where you want figure to be saved
% whichCombi: 1 = compare to "All" maps; 2 = compare to most overlapping
%       map; 3 = compare to same map; 4 = compare to all overlapping maps;
%       5 = use voxel selection from other data structure (both use whichCombi = 4) 
%
% Output
% saves info in stat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minAmountDataPoints = 1;
n_fig = 0;

%% load relevant statistics
if whichCombi == 5
    num_stat_numEven_timeOdd = load('stat_NumerosityEven_NumerosityOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=5_selection=NumerosityEven_TimingOdd_topo_measure=mean.mat');
    num_stat_numOdd_timeEven = load('stat_NumerosityEven_NumerosityOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=5_selection=NumerosityOdd_TimingEven_topo_measure=mean.mat');
    time_stat_numEven_timeOdd = load('stat_TimingEven_TimingOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=5_selection=NumerosityEven_TimingOdd_topo_measure=mean.mat');
    time_stat_numOdd_timeEven = load('stat_TimingEven_TimingOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=5_selection=NumerosityOdd_TimingEven_topo_measure=mean.mat');

elseif whichCombi == 4
    num_stat = load('stat_NumerosityEven_NumerosityOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat');
    time_stat = load('stat_TimingEven_TimingOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat');
    stat_numEven_timeOdd = load('stat_NumerosityEven_TimingOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat');
    stat_numOdd_timeEven = load('stat_NumerosityOdd_TimingEven_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat');
else
    return
end
clear stat

% Num-Dur and Num-Per
parameters_diff = {'x0.x0', 'x0.y0'};

% NumEven and TimeEven should be compared to NumOdd-NumEven, but
% voxel selected based on the Even runs

%% calculate mean correlations 
if whichCombi == 5
    num_corrs_x0x0 = NumTim_make_mean_correlation(num_stat_numEven_timeOdd.stat.data.NumerosityEven.correlation.x0.x0,num_stat_numOdd_timeEven.stat.data.NumerosityOdd.correlation.x0.x0);
    time_corrs_x0x0 = NumTim_make_mean_correlation(time_stat_numEven_timeOdd.stat.data.TimingEven.correlation.x0.x0,time_stat_numOdd_timeEven.stat.data.TimingOdd.correlation.x0.x0);
    time_corrs_y0y0 = NumTim_make_mean_correlation(time_stat_numEven_timeOdd.stat.data.TimingEven.correlation.y0.y0,time_stat_numOdd_timeEven.stat.data.TimingOdd.correlation.y0.y0);
    num_maps = num_stat_numOdd_timeEven.stat.data.NumerosityEven.map;
    time_maps = time_stat_numOdd_timeEven.stat.data.TimingEven.map;

    diff_corrs_x0x0_num = NumTim_make_mean_correlation(num_stat_numEven_timeOdd.stat.data.NumerosityEven.selection_correlation.x0.x0,num_stat_numOdd_timeEven.stat.data.NumerosityOdd.selection_correlation.x0.x0);
    diff_corrs_x0y0_num = NumTim_make_mean_correlation(num_stat_numEven_timeOdd.stat.data.NumerosityEven.selection_correlation.x0.y0,num_stat_numOdd_timeEven.stat.data.NumerosityOdd.selection_correlation.x0.y0);
    diff_corrs_x0x0_time = NumTim_make_mean_correlation(time_stat_numEven_timeOdd.stat.data.TimingOdd.selection_correlation.x0.x0,time_stat_numOdd_timeEven.stat.data.TimingEven.selection_correlation.x0.x0);
    diff_corrs_x0y0_time = NumTim_make_mean_correlation(time_stat_numEven_timeOdd.stat.data.TimingOdd.selection_correlation.x0.y0,time_stat_numOdd_timeEven.stat.data.TimingEven.selection_correlation.x0.y0);

    diff_maps_num = num_stat_numEven_timeOdd.stat.data.NumerosityEven.map;
    diff_maps_time = time_stat_numEven_timeOdd.stat.data.TimingOdd.map;

    subj_num = num_stat_numEven_timeOdd.stat.data.NumerosityEven.subj;
    hemi_num = num_stat_numEven_timeOdd.stat.data.NumerosityEven.hemi;

    subj_time = time_stat_numEven_timeOdd.stat.data.TimingOdd.subj;
    hemi_time = time_stat_numEven_timeOdd.stat.data.TimingOdd.hemi;

elseif whichCombi == 4
    num_corrs_x0x0 = NumTim_make_mean_correlation(num_stat.stat.data.NumerosityEven.correlation.x0.x0,num_stat.stat.data.NumerosityOdd.correlation.x0.x0);
    time_corrs_x0x0 = NumTim_make_mean_correlation(time_stat.stat.data.TimingEven.correlation.x0.x0,time_stat.stat.data.TimingOdd.correlation.x0.x0);
    time_corrs_y0y0 = NumTim_make_mean_correlation(time_stat.stat.data.TimingEven.correlation.y0.y0,time_stat.stat.data.TimingOdd.correlation.y0.y0);
    num_maps = num_stat.stat.data.NumerosityEven.map;
    time_maps = time_stat.stat.data.TimingEven.map;


    diff_stat_1 = stat_numEven_timeOdd.stat;
    diff_stat_2 = stat_numOdd_timeEven.stat;

    diff_corrs_x0x0_num = NumTim_make_mean_correlation(diff_stat_1.data.NumerosityEven.correlation.x0.x0,diff_stat_2.data.NumerosityOdd.correlation.x0.x0);
    diff_corrs_x0y0_num = NumTim_make_mean_correlation(diff_stat_1.data.NumerosityEven.correlation.x0.y0,diff_stat_2.data.NumerosityOdd.correlation.x0.y0);
    diff_corrs_x0x0_time = NumTim_make_mean_correlation(diff_stat_1.data.TimingOdd.correlation.x0.x0,diff_stat_2.data.TimingEven.correlation.x0.x0);
    diff_corrs_x0y0_time = NumTim_make_mean_correlation(diff_stat_1.data.TimingOdd.correlation.x0.y0,diff_stat_2.data.TimingEven.correlation.x0.y0);

    diff_maps_num = diff_stat_1.data.NumerosityEven.map;
    diff_maps_time = diff_stat_1.data.TimingOdd.map;

    subj_num = diff_stat_1.data.NumerosityEven.subj;
    hemi_num = diff_stat_1.data.NumerosityEven.hemi;

    subj_time = diff_stat_1.data.TimingOdd.subj;
    hemi_time = diff_stat_1.data.TimingOdd.hemi;

end
time_rois = ["TLO","TTOP","TTOA","TPO","TLS","TPCI","TPCM","TPCS","TFI","TFS"];
num_rois = ["NLO","NTO","NPO","NPCI","NPCM","NPCS","NFI","NFS"];


%% check normality
DT_runs = {'num', 'time'};
diff_parameters = {'x0x0','x0y0'};
num_parameters = {'x0x0', 'x0x0'};
time_parameters = {'x0x0','y0y0'};
p_norms = [];

for run = 1:length(DT_runs)
    % manually checked: same order everywhere
    eval(['diff_maps = diff_maps_', DT_runs{run},';'])

    for par = 1:length(diff_parameters)
        eval(['diff_corrs = diff_corrs_', diff_parameters{par}, '_', DT_runs{run},';'])
        eval(['same_parameters =', DT_runs{run}, '_parameters;'])

        eval(['same_corrs = ',DT_runs{run},'_corrs_', same_parameters{par},';'])
        eval(['rois =', DT_runs{run},'_rois;'])

        for roi = 1:length(rois)
            if sum(~isnan(same_corrs(diff_maps==rois{roi}))) > minAmountDataPoints
                try
                    [~,stat.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).(rois{roi}).normality,~] = force_swtest(same_corrs(diff_maps==rois{roi}) - diff_corrs(diff_maps==rois{roi}),0.05,1);
                catch
                    stat.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).(rois{roi}).normality = 0;
                end
                p_norms = [p_norms stat.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).(rois{roi}).normality];

            end
        end

    end
end
disp(p_norms)

%% pairwise comparisons
close all
nboot = 1000;


for run = 1:length(DT_runs)
    eval(['diff_maps = diff_maps_', DT_runs{run},';'])


    % FDR correction takes place within numerosity maps when
    % NumEven-DurEven vs. NumEven-NumOdd and NumEven-PerEven vs. NumEven-NumOdd
    % NOT INCLUDED THERE:
    % NumOdd-DurOdd and NumOdd-PerOdd, since these are done on
    % different voxel selections (overlap with NumEven-TimeEven and
    % NumOdd-TimeOdd, respectively)
    pvals= [];
    store_rois = [];
    comparison = [];

    for par = 1:length(diff_parameters)
        eval(['diff_corrs = diff_corrs_', diff_parameters{par}, '_', DT_runs{run},';'])
        eval(['same_parameters =', DT_runs{run}, '_parameters;'])

        eval(['same_corrs = ',DT_runs{run},'_corrs_', same_parameters{par},';'])
        eval(['rois =', DT_runs{run},'_rois;'])

        % so you can make paired comparisons
        % already the case for whichCombi == 5
        % just for plots, wilcoxon does it right automatically
        same_corrs(isnan(diff_corrs)) = NaN;
        diff_corrs(isnan(same_corrs)) = NaN;

        x_counter = 0;
        pos_rois = [];
        figure();
        for roi = 1:length(rois)
            if sum(~isnan(same_corrs(diff_maps==rois{roi}) - diff_corrs(diff_maps==rois{roi}))) > minAmountDataPoints
                [stat.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).(rois{roi}).p,stat.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).(rois{roi}).h,stat.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).(rois{roi}).stats] = signrank(same_corrs(diff_maps==rois{roi}),diff_corrs(diff_maps==rois{roi}),'method','approximate');%,'tail','right');

                % summary stats with their CI
                x_counter = x_counter + 1;
                point_x_same(x_counter) = x_counter + 0.2-0.5;
                point_x_diff(x_counter) = x_counter - 0.2-0.5;

                point_y_same(x_counter) = nanmedian(same_corrs(diff_maps==rois{roi}));
                point_y_diff(x_counter) = nanmedian(diff_corrs(diff_maps==rois{roi}));

                bootstr_CI_same = bootci(nboot, {@nanmedian, same_corrs(diff_maps==rois{roi})},'alpha',0.05);
                error_bars_same(:,x_counter) = [point_y_same(x_counter)-bootstr_CI_same(1,:), bootstr_CI_same(2,:)-point_y_same(x_counter)];

                bootstr_CI_diff = bootci(nboot, {@nanmedian, diff_corrs(diff_maps==rois{roi})},'alpha',0.05);
                error_bars_diff(:,x_counter) = [point_y_diff(x_counter)-bootstr_CI_diff(1,:), bootstr_CI_diff(2,:)-point_y_diff(x_counter)];

                pos_rois = [pos_rois rois(roi)];

            else
                stat.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).(rois{roi}).p = NaN;

            end

            pvals = [pvals stat.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).(rois{roi}).p];
            comparison = [comparison same_parameters(par)];
            store_rois = [store_rois rois(roi)];

        end

        size_figure = 9; % 90 single column, 140 1.5, 190 double
        n_fig = n_fig + 1;
        figure(n_fig);
        hold on

        % plot information
        hold on
        colorMapScatter = linspace(1,10,8);
        cm = getPyPlot_cMap('Set2');
        colormap(gca, cm)


        eval(['subj = subj_',DT_runs{run},';'])
        eval(['hemi = hemi_',DT_runs{run},';'])

        subjs = str2num(char(strip(subj,'S')));
        [~, pos]= ismember(diff_maps,pos_rois);
        % these data points are not compared, so should not be
        % in the plot
        same_corrs(pos==0)=NaN;
        diff_corrs(pos==0)=NaN;

        left_color = colorMapScatter(subjs(hemi=='Left'));
        right_color = colorMapScatter(subjs(hemi=='Right'));


        fig = gca;
        fig.Units = "centimeters";
        fig.Color = 'w';
        set(gca, 'Fontsize',8,'TickDir','out', 'FontName','Arial','Units','centimeters','Color','w','XColor','k','YColor','k','box','off');
        fig.YAxis.TickLabelColor = [0 0 0];
        fig.YAxis.FontName = 'Arial';
        fig.YAxis.FontSize = 7;
        fig.YAxis.Color = [0 0 0];
        fig.XAxis.TickLabelColor = [0 0 0];
        fig.XAxis.FontName = 'Arial';
        fig.XAxis.FontSize = 7;
        fig.XAxis.Color = [0 0 0];
        fig.XTickLabelRotation = 90;
        size_axis1 = (size_figure-2)/(8+ 10)*8;
        size_axis2 = (size_figure-2)/(8+ 10)*10;
        height_figs = size_axis2;

        size_axis = eval(['size_axis',char(num2str(run)),';']);
        fig.Position = [1 1 size_axis height_figs];
        fig.YLabel.String = 'Correlation preferred numerosity and time';
        fig.XLabel.String = [DT_runs{run}, ' maps'];


        scatter(pos(hemi=='Left')+0.2 -0.5 -0.05,same_corrs(hemi=='Left'), 2.5*1.4,left_color(:),'o','filled');
        scatter(pos(hemi=='Right')+ 0.2 -0.5 +0.05, same_corrs(hemi=='Right'), 2.5*1.4,right_color(:),'d','filled');
        scatter(pos(hemi=='Left')- 0.2 -0.5 -0.05,diff_corrs(hemi=='Left'), 2.5*1.4,left_color(:),'o','filled');
        scatter(pos(hemi=='Right')- 0.2 -0.5+0.05, diff_corrs(hemi=='Right'),2.5*1.4,right_color(:),'d','filled');
        ylim = [-1,1];
        plot([0 length(pos_rois)+size_axis/length(pos_rois)*0.5],[0 0], 'k-')



        e = errorbar(point_x_same,point_y_same,error_bars_same(1,:),error_bars_same(2,:),'s','CapSize',0,'LineWidth',1,'MarkerFaceColor','k','MarkerEdgeColor','None','MarkerSize',3);
        e.Color = [0 0 0];
        e = errorbar(point_x_diff,point_y_diff,error_bars_diff(1,:),error_bars_diff(2,:),'s','CapSize',0,'LineWidth',1,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','None','MarkerSize',3);
        e.Color = [0.5 0.5 0.5];

        ax = gca;
        set(ax,'XLim',[0,length(pos_rois)+size_axis/length(pos_rois)*0.5],'YLim',ylim,'XTick',0.5:length(pos_rois)-0.5,'XTickLabel',pos_rois)

        plot_values.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).errorbars.same = error_bars_same;
        plot_values.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).errorbars.diff = error_bars_diff;
        plot_values.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).point.same = point_y_same;
        plot_values.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).point.diff = point_y_diff;
        plot_values.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).x.same = point_x_same;
        plot_values.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).x.diff = point_x_diff;
        plot_values.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).roi = pos_rois;
        plot_values.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).scatter.corr.same = same_corrs;
        plot_values.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).scatter.corr.diff = diff_corrs;
        plot_values.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).scatter.maps = diff_maps;
        plot_values.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).scatter.subj = subj;
        plot_values.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).scatter.hemi = hemi;

        savename_eps = strcat(fig_path, '_',DT_runs{run}, '_',diff_parameters{par},'_',same_parameters{par},'.eps');
        disp(savename_eps)
        %         saveas(gcf,savename_eps,'epsc')
        close all
    end

    % check whether they go here
    adj_p = NaN(1,length(pvals));
    [~,~,~,adj_p(~isnan(pvals))] = fdr_bh(pvals(~isnan(pvals)));

    for par = 1:length(diff_parameters)
        eval(['same_parameters =', DT_runs{run}, '_parameters;'])
        for roi = 1:length(rois)
            stat.(DT_runs{run}).correlation.(diff_parameters{par}).(same_parameters{par}).(rois{roi}).adj_p = adj_p(roi+length(rois)*(par-1));
        end
    end
    disp(adj_p)
    disp(DT_runs{run})
    disp(store_rois(adj_p>=0.05))
    disp(comparison(adj_p>=0.05))
    disp(adj_p(adj_p>=0.05))

end


savename_struct = strcat(fig_path, '_plot_values.mat');
save(savename_struct, 'plot_values')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% statistics against chance  FOR REPEATED MEASURES

p_norms = [];
for run = 1:length(DT_runs)
    eval(['diff_maps = diff_maps_', DT_runs{run},';'])

    store_rois = [];
    comparison = [];
    % done for duration AND period in same maps in order to apply FDR
    pvals_chance = [];

    eval(['rois =', DT_runs{run},'_rois;'])
    eval(['same_parameters =', DT_runs{run}, '_parameters;'])



    for par = 1:length(same_parameters)

        % This is for repeated measures
        % For different num-time runs voxel selection will stay the same
        % for numerosity num-num is used for both
        % for period and duration comparisons
        % so only 8 different values
        if DT_runs{run} == "num" && par == 2
            continue
        end

        eval(['diff_corrs = diff_corrs_', diff_parameters{par}, '_', DT_runs{run},';'])
        eval(['same_corrs = ',DT_runs{run},'_corrs_', same_parameters{par},';'])


        %% Dont (officially) need to do this for chance calculation, but now it includes the same datapoints as the paired comparison
        % this makes a difference in whichCombi == 5
        same_corrs(isnan(diff_corrs)) = NaN;
        diff_corrs(isnan(same_corrs)) = NaN;

        for roi = 1:length(rois)
            try
                [~,stat.(DT_runs{run}).chance_distribution.correlation.(same_parameters{par}).(rois{roi}).normality,~] = force_swtest(same_corrs(diff_maps==rois{roi}),0.05,1);
            catch
                stat.(DT_runs{run}).chance_distribution.correlation.(same_parameters{par}).(rois{roi}).normality = 0;
            end
            p_norms = [p_norms stat.(DT_runs{run}).chance_distribution.correlation.(same_parameters{par}).(rois{roi}).normality];


            if sum(~isnan(same_corrs(diff_maps==rois{roi}))) > minAmountDataPoints
                [stat.(DT_runs{run}).chance_distribution.correlation.(same_parameters{par}).(rois{roi}).p,stat.(DT_runs{run}).chance_distribution.correlation.(same_parameters{par}).(rois{roi}).h,stat.(DT_runs{run}).chance_distribution.correlation.(same_parameters{par}).(rois{roi}).stats] = signrank(same_corrs(diff_maps==rois{roi}),0,'method','approximate');%,'tail','right');

            else
                stat.(DT_runs{run}).chance_distribution.correlation.(same_parameters{par}).(rois{roi}).p = NaN;
            end
            pvals_chance = [pvals_chance stat.(DT_runs{run}).chance_distribution.correlation.(same_parameters{par}).(rois{roi}).p];
            store_rois = [store_rois rois(roi)];
            comparison = [comparison same_parameters(par)];
        end
    end

    adj_p = nan(size(pvals_chance));
    [~,~,~,adj_p(~isnan(pvals_chance))] = fdr_bh(pvals_chance(~isnan(pvals_chance)));

    for par = 1:length(same_parameters)
        if DT_runs{run} == "num" && par == 2
            for roi = 1:length(rois)
                stat.(DT_runs{run}).chance_distribution.correlation.(same_parameters{par}).(rois{roi}).adj_p = adj_p(roi);
            end
        else
            for roi = 1:length(rois)
                stat.(DT_runs{run}).chance_distribution.correlation.(same_parameters{par}).(rois{roi}).adj_p = adj_p(roi+length(rois)*(par-1));
            end
        end
    end

end


savename_struct = strcat(fig_path, '_stat.mat');
save(savename_struct, 'stat')