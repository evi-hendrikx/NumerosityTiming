function NumTim_compare_proportions(fig_path)
%% Compares proportions for Num-Time vs repeated measures
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

% will not compare if there's only one data point in both conditions
minAmountDataPoints = 1;
n_fig = 0;

% load in necessary data
num_stat = load('stat_NumerosityEven_NumerosityOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat');
time_stat = load('stat_TimingEven_TimingOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat');
stat_numEven_timeOdd = load('stat_NumerosityEven_TimingOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat');
stat_numOdd_timeEven = load('stat_NumerosityOdd_TimingEven_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat');
clear stat

% NumEven and TimeEven should be compared to NumOdd-NumEven, but
% voxel selected based on the Even runs
num_props = nanmean([num_stat.stat.data.NumerosityEven.prop,num_stat.stat.data.NumerosityOdd.prop],2);
time_props = nanmean([time_stat.stat.data.TimingEven.prop,time_stat.stat.data.TimingOdd.prop],2);
diff_props_num = nanmean([stat_numEven_timeOdd.stat.data.NumerosityEven.prop,stat_numOdd_timeEven.stat.data.NumerosityOdd.prop],2);
diff_props_time = nanmean([stat_numEven_timeOdd.stat.data.TimingOdd.prop,stat_numOdd_timeEven.stat.data.TimingEven.prop],2);

num_maps = num_stat.stat.data.NumerosityEven.map;
time_maps = time_stat.stat.data.TimingEven.map;
subj_num = stat_numEven_timeOdd.stat.data.NumerosityEven.subj;
hemi_num = stat_numEven_timeOdd.stat.data.NumerosityEven.hemi;
subj_time = stat_numEven_timeOdd.stat.data.TimingOdd.subj;
hemi_time = stat_numEven_timeOdd.stat.data.TimingOdd.hemi;

time_rois = ["TLO","TTOP","TTOA","TPO","TLS","TPCI","TPCM","TPCS","TFI","TFS"];
num_rois = ["NLO","NTO","NPO","NPCI","NPCM","NPCS","NFI","NFS"];

%% check normality
DT_runs = {'num', 'time'};
GNames = {num_rois,time_rois};
p_norms = [];
for run = 1:length(DT_runs)
    % manually checked: same order everywhere
    eval(['diff_maps =', DT_runs{run}, '_maps;'])
    eval(['diff_props = diff_props_', DT_runs{run},';'])
    eval(['same_props = ',DT_runs{run},'_props;'])
    eval(['rois =', DT_runs{run},'_rois;'])

    for roi = 1:length(rois)
        % manually check TODO
        if sum(~isnan(same_props(diff_maps==rois{roi}))) > minAmountDataPoints
            try
                [~,stat.(DT_runs{run}).proportion.(rois{roi}).normality,~] = force_swtest(same_props(diff_maps==rois{roi}) - diff_props(diff_maps==rois{roi}),0.05,1);
            catch
                stat.(DT_runs{run}).proportion.(rois{roi}).normality = 0;
            end
            p_norms = [p_norms stat.(DT_runs{run}).proportion.(rois{roi}).normality];

        end
    end

end
disp(p_norms)

%% check pairwise comparisons
close all
nboot = 1000;


for run = 1:length(DT_runs)
    % FDR correction takes place within numerosity maps and timing maps
    % separately
    pvals= [];
    store_rois = [];
    comparison = [];

    eval(['diff_maps = ' DT_runs{run},'_maps;'])
    eval(['diff_props = diff_props_', DT_runs{run},';'])
    eval(['same_props = ',DT_runs{run},'_props;'])
    eval(['rois =', DT_runs{run},'_rois;'])

    % so you can make paired comparisons
    same_props(isnan(diff_props)) = NaN;
    diff_props(isnan(same_props)) = NaN;

    x_counter = 0;
    pos_rois = [];
    figure();

    for roi = 1:length(rois)

        if sum(~isnan(same_props(diff_maps==rois{roi}) - diff_props(diff_maps==rois{roi}))) > minAmountDataPoints
            [stat.(DT_runs{run}).proportion.(rois{roi}).p,stat.(DT_runs{run}).proportion.(rois{roi}).h,stat.(DT_runs{run}).proportion.(rois{roi}).stats] = signrank(same_props(diff_maps==rois{roi}),diff_props(diff_maps==rois{roi}),'method','approximate');%,'tail','right');

            % summary stats with their CI
            x_counter = x_counter + 1;
            point_x_same(x_counter) = x_counter + 0.2 -0.5;
            point_x_diff(x_counter) = x_counter - 0.2 -0.5;

            point_y_same(x_counter) = nanmedian(same_props(diff_maps==rois{roi}));
            point_y_diff(x_counter) = nanmedian(diff_props(diff_maps==rois{roi}));

            bootstr_CI_same = bootci(nboot, {@nanmedian, same_props(diff_maps==rois{roi})},'alpha',0.05);
            error_bars_same(:,x_counter) = [point_y_same(x_counter)-bootstr_CI_same(1,:), bootstr_CI_same(2,:)-point_y_same(x_counter)];

            bootstr_CI_diff = bootci(nboot, {@nanmedian, diff_props(diff_maps==rois{roi})},'alpha',0.05);
            error_bars_diff(:,x_counter) = [point_y_diff(x_counter)-bootstr_CI_diff(1,:), bootstr_CI_diff(2,:)-point_y_diff(x_counter)];

            pos_rois = [pos_rois rois(roi)];

        else
            stat.(DT_runs{run}).proportion.(rois{roi}).p = NaN;

        end

        pvals = [pvals stat.(DT_runs{run}).proportion.(rois{roi}).p];
        store_rois = [store_rois rois(roi)];

    end
    size_figure = 9; % 90 single column, 140 1.5, 190 double
    n_fig = n_fig + 1;
    figure(n_fig);
    hold on


    % plot information
    colorMapScatter = linspace(1,10,8);
    cm = getPyPlot_cMap('Set2');
    colormap(gca, cm);

    eval(['subj = subj_',DT_runs{run},';']);
    eval(['hemi = hemi_',DT_runs{run},';']);


    subjs = str2num(char(strip(subj,'S')));
    [~, pos]= ismember(diff_maps,pos_rois);

    % these data points are not compared, so should not be
    % in the plot
    same_props(pos==0)=NaN;
    diff_props(pos==0)=NaN;

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
    size_axis1 = (size_figure-2)/(length(GNames{1})+ length(GNames{2}))*length(GNames{1});
    size_axis2 = (size_figure-2)/(length(GNames{1})+ length(GNames{2}))*length(GNames{2});
    if length(GNames{1}) > length(GNames{2})
        height_figs = size_axis1;
    else
        height_figs = size_axis2;
    end
    size_axis = eval(['size_axis',char(num2str(run)),';']);
    fig.Position = [1 1 size_axis height_figs];
    fig.YLabel.String = 'Proportion of overlap';
    fig.XLabel.String = [DT_runs{run}, ' maps'];


    scatter(pos(hemi=='Left')+ 0.2 -0.5 -0.05,same_props(hemi=='Left'), 2.5*1.4,left_color(:),'o','filled');
    scatter(pos(hemi=='Right')+ 0.2 -0.5 +0.05, same_props(hemi=='Right'), 2.5*1.4,right_color(:),'d','filled');
    scatter(pos(hemi=='Left')- 0.2 -0.5 -0.05,diff_props(hemi=='Left'), 2.5*1.4,left_color(:),'o','filled');
    scatter(pos(hemi=='Right')- 0.2 -0.5+0.05, diff_props(hemi=='Right'),2.5*1.4,right_color(:),'d','filled');


    e = errorbar(point_x_same,point_y_same,error_bars_same(1,:),error_bars_same(2,:),'s','CapSize',0,'LineWidth',1,'MarkerFaceColor','k','MarkerEdgeColor','None','MarkerSize',3);
    e.Color = [0 0 0];
    e = errorbar(point_x_diff,point_y_diff,error_bars_diff(1,:),error_bars_diff(2,:),'s','CapSize',0,'LineWidth',1,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','None','MarkerSize',3);
    e.Color = [0.5 0.5 0.5];

    ax = gca;
    set(ax,'XLim',[0,length(GNames{run})+size_axis/length(GNames{run})*0.5],'YLim',ylim,'XTick',0.5:length(GNames{run})-0.5,'XTickLabel',GNames{run})

    plot_values.(DT_runs{run}).proportion.errorbars.same = error_bars_same;
    plot_values.(DT_runs{run}).proportion.errorbars.diff = error_bars_diff;
    plot_values.(DT_runs{run}).proportion.point.same = point_y_same;
    plot_values.(DT_runs{run}).proportion.point.diff = point_y_diff;
    plot_values.(DT_runs{run}).proportion.x.same = point_x_same;
    plot_values.(DT_runs{run}).proportion.x.diff = point_x_diff;
    plot_values.(DT_runs{run}).proportion.roi = pos_rois;
    plot_values.(DT_runs{run}).proportion.scatter.prop.same = same_props;
    plot_values.(DT_runs{run}).proportion.scatter.prop.diff = diff_props;
    plot_values.(DT_runs{run}).proportion.scatter.maps = diff_maps;
    plot_values.(DT_runs{run}).proportion.scatter.subj = subj;
    plot_values.(DT_runs{run}).proportion.scatter.hemi = hemi;

    % check whether they go here
    adj_p = NaN(1,length(pvals));
    [~,~,~,adj_p(~isnan(pvals))] = fdr_bh(pvals(~isnan(pvals)));
    for roi = 1:length(rois)
        stat.(DT_runs{run}).proportion.(rois{roi}).adj_p = adj_p(roi);
    end
    disp(adj_p)
    disp(DT_runs{run})
    disp(store_rois(adj_p>=0.05))
    disp(adj_p(adj_p>=0.05))

    savename_eps = strcat(fig_path, '_',DT_runs{run}, '.eps');
    disp(savename_eps)
    saveas(gcf,savename_eps,'epsc')
    close all

end

savename_struct = strcat(fig_path, '_plot_values.mat');
save(savename_struct, 'plot_values')

savename_struct = strcat(fig_path, '_stat.mat');
save(savename_struct, 'stat')

