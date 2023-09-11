function NumTim_plots_prop_corr(fig_path, NumTim_data, stat, whichCombi, DT_runs)
%% Make boxplot of the requested results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% fig_path: where you want to save the figure
% NumTim_data: made with the NumTim_load_data script. Contains all values
%       that you could need for all following analyses
% stat: structure containing statistical info from all the performed tests
% whichCombi: 1 = compare to "All" maps; 2 = compare to most overlapping
%       map; 3 = compare to same map; 4 = all overlapping maps
% DT_runs: the two data runs you want to compare (options: TimingAll, TimingEven,
%       TimingOdd, NumerosityAll, NumerosityEven, NumerosityOdd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
n_fig = 0;
minAmountDataPoints = 1;

%% get information
[~,~,Run_1, Run_2] = NumTim_get_info(NumTim_data,DT_runs);
GNames = {{Run_1.Names{1:end-1}},{Run_2.Names{1:end-1}}};

variable = {'proportion','correlation'};
for var = 1:length(variable)
    for run = 1:length(DT_runs)

        % prepare information
        parameters = []; data = [];
        if variable{var} == "proportion"
            data = {stat.data.(DT_runs{run}).prop};
            if whichCombi == 1 || whichCombi == 4
                chance_data = {stat.data.(DT_runs{run}).chance_proportion};
            end
        elseif variable{var} == "correlation"
            parameters = [parameters {'x0x0'}];
            data = {stat.data.(DT_runs{run}).correlation.x0.x0};
            if isfield(stat.data.(DT_runs{run}).(variable{var}), 'y0') && string(Run_1.Type) ~= string(Run_2.Type) % only want to know dur-dur per-per and num-num if even vs odd runs
                parameters = [parameters {'y0x0'}];
                data = [data {stat.data.(DT_runs{run}).correlation.y0.x0}];
            end
            if isfield(stat.data.(DT_runs{run}).(variable{var}).x0, 'y0') && string(Run_1.Type) ~= string(Run_2.Type) % only want to know dur-dur per-per and num-num if even vs odd runs
                parameters = [parameters {'x0y0'}];
                data = [data {stat.data.(DT_runs{run}).correlation.x0.y0}];
            end
            if isfield(stat.data.(DT_runs{run}).(variable{var}), 'y0') && isfield(stat.data.(DT_runs{run}).(variable{var}).x0, 'y0')
                parameters = [parameters {'y0y0'}];
                data = [data {stat.data.(DT_runs{run}).correlation.y0.y0}];
            end
        end

        if variable{var} =="correlation"
            %     % kick out areas with only 1 datapoint from correlational
            %     % analyses
            loop_roi_names = GNames{run};
            kicked_rois = 0;
            for roi = 1:length(loop_roi_names)
                for par = 1:length(parameters)
                    if eval(['sum(~isnan(data{par}(stat.data.(DT_runs{run}).map == loop_roi_names{roi}))) <= minAmountDataPoints'])
                        data{par}(stat.data.(DT_runs{run}).map == loop_roi_names{roi}) = NaN;
%                         if par == 1
%                             GNames{run}(roi - kicked_rois) = [];
%                             kicked_rois = kicked_rois + 1;
%                         end
                    end
                end
            end

        end
        

        %% properties of figures
        size_figure = 9; % 90 single column, 140 1.5, 190 double

        for par = 1:length(data)
            clear bootstap error_bars
            n_fig = n_fig + 1;
            figure(n_fig);
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

            %% determine properties subplots
            if string(variable{var}) == "proportion"
                fig.YLabel.String = 'Proportion of overlap';
                if string(Run_1.Type) ~= string(Run_2.Type)
                    eval(['fig.XLabel.String = [Run_',char(num2str(run)),'.Type ', char("' maps'"),'];']);
                else
                    eval(['fig.XLabel.String = [Run_',char(num2str(run)),'.Type Run_', char(num2str(run)),'.Run ', char("' runs'"),']']);
                end

            elseif string(Run_1.Type) == string(Run_2.Type)  && string(variable{var}) == "correlation" % only dur-dur and per-per are necessary
                eval(['fig.XLabel.String = [Run_',char(num2str(run)),'.Type ', char("' maps'"),'];']);

                if string(Run_1.Type) == "Timing"
                    if par == 1
                        fig.YLabel.String = ['Correlation preferred duration (rs)'];
                    else
                        fig.YLabel.String = ['Correlation preferred period (rs)'];
                    end
                else
                    fig.YLabel.String = 'Correlation preferred numerosity (rs)';
                end

            elseif string(variable{var}) == "correlation"
                eval(['fig.XLabel.String = [Run_',char(num2str(run)),'.Type ', char("' maps'"),'];']);

                if par == 1
                    fig.YLabel.String = ['Correlation preferred numerosity and duration (rs)'];
                else
                    fig.YLabel.String = ['Correlation preferred numerosity and period (rs)'];
                end
            end

            point_x=[];point_y=[];z_mean=[];point_y_median=[];


            % summary stats with their CI
            nboot = 1000;
            for roi = 1:length(GNames{run})


                point_x(roi) = roi-0.5;
                if (string(variable(var)) == "proportion" && (stat.ANOVA.(DT_runs{1}).(variable{var}).p_norm >= 0.05 && stat.ANOVA.(DT_runs{2}).(variable{var}).p_norm >= 0.05)) || ...
                        (string(variable{var}) == "correlation" && (stat.ANOVA.(DT_runs{1}).(variable{var}).x0x0.p_norm >= 0.05 && stat.ANOVA.(DT_runs{2}).(variable{var}).x0x0.p_norm >= 0.05) && ...
                        whichCombi ~= 5 && ((string(Run_1.Type) ~= string(Run_2.Type) && stat.ANOVA.(DT_runs{1}).(variable{var}).x0y0.p_norm >= 0.05 && stat.ANOVA.(DT_runs{2}).(variable{var}).x0y0.p_norm >= 0.05)...
                        ||(string(Run_1.Type) == string(Run_2.Type) && string(Run_2.Type) =="Timing" && stat.ANOVA.(DT_runs{1}).(variable{var}).y0y0.p_norm >= 0.05)))
                    point_y(roi) = nanmean(data{par}(stat.data.(DT_runs{run}).map == GNames{run}(roi)));
                    if length(data{par}(stat.data.(DT_runs{run}).map == GNames{run}(roi))) > 1
                        SEM = nanstd(data{par}(stat.data.(DT_runs{run}).map == GNames{run}(roi)))/ sqrt(sum(~isnan(data{par}(stat.data.(DT_runs{run}).map == GNames{run}(roi)))));
                        error_bars(:,roi) = [SEM SEM]';
                    else
                        error_bars(:,roi) = [NaN NaN];
                    end
                else

                    point_y(roi) = nanmedian(data{par}(stat.data.(DT_runs{run}).map == GNames{run}(roi)));
                    if sum(~isnan(data{par}(stat.data.(DT_runs{run}).map == GNames{run}(roi)))) > minAmountDataPoints
                        bootstr_CI = bootci(nboot, {@nanmedian, data{par}(stat.data.(DT_runs{run}).map == GNames{run}(roi))},'alpha',0.05);
                        error_bars(:,roi) = [point_y(roi)-bootstr_CI(1,:), bootstr_CI(2,:)-point_y(roi)];
                    else
                        error_bars(:,roi) = [NaN NaN];
                    end

                end
            end


            hold on
            colorMapScatter = linspace(1,10,8);
            cm = getPyPlot_cMap('Set2');
            colormap(gca, cm)
            subjs = str2num(char(strip(stat.data.(DT_runs{run}).subj,'S')));
            left_color = colorMapScatter(subjs(stat.data.(DT_runs{run}).hemi=='Left'));
            right_color = colorMapScatter(subjs(stat.data.(DT_runs{run}).hemi=='Right'));


            if string(variable{var}) == "proportion" || string(variable{var}) == "correlation"


                [~, pos] = ismember(stat.data.(DT_runs{run}).map,GNames{run});
                pos = pos - 0.5;


                scatter(pos(stat.data.(DT_runs{run}).hemi=='Left')-0.15,data{par}(stat.data.(DT_runs{run}).hemi=='Left'), 2.5*1.4,left_color(:),'o','filled');
                scatter(pos(stat.data.(DT_runs{run}).hemi=='Right')+0.15, data{par}(stat.data.(DT_runs{run}).hemi=='Right'), 2.5*1.4,right_color(:),'d','filled');

                if string(variable{var}) == "proportion"
                    ylim = [0,1];
                    if  whichCombi == 1||whichCombi == 4
                        scatter(pos(stat.data.(DT_runs{run}).hemi=='Left')-0.15, chance_data{par}(stat.data.(DT_runs{run}).hemi=='Left'),1.3*1.4,left_color(:),'_','LineWidth',height_figs/15);
                        scatter(pos(stat.data.(DT_runs{run}).hemi=='Right')+0.15, chance_data{par}(stat.data.(DT_runs{run}).hemi=='Right'),1.3*1.4,right_color(:),'_','LineWidth',height_figs/15);
                    end
                elseif string(variable{var}) == "correlation"
                    ylim = [-1,1];
                    plot([0 length(GNames{run})-0.5+size_axis/length(GNames{run})*0.5],[0 0], 'k-')
                end

                e = errorbar(point_x,point_y,error_bars(1,:),error_bars(2,:),'s','CapSize',0,'LineWidth',1,'MarkerFaceColor','k','MarkerEdgeColor','None','MarkerSize',3);
                e.Color = [0 0 0];

                ax = gca;


                set(ax,'XLim',[0,length(GNames{run})-0.5+size_axis/length(GNames{run})*0.5],'YLim',ylim,'XTick',0.5:length(GNames{run})-0.5,'XTickLabel',GNames{run})

                if string(variable{var}) == "correlation"
                    plot_values.(DT_runs{run}).(variable{var}).(parameters{par}).errorbars = error_bars;
                    plot_values.(DT_runs{run}).(variable{var}).(parameters{par}).point = point_y;
                    plot_values.(DT_runs{run}).(variable{var}).(parameters{par}).x = point_x;
                elseif string(variable{var}) == "proportion"

                    plot_values.(DT_runs{run}).(variable{var}).errorbars = error_bars;
                    plot_values.(DT_runs{run}).(variable{var}).point = point_y;
                    plot_values.(DT_runs{run}).(variable{var}).x = point_x;
                end
           
            end
            try
                savename_eps = strcat(fig_path, '_', variable{var}, '_',DT_runs{run},'_',parameters{par}, '.eps');
            catch
                savename_eps = strcat(fig_path, '_', variable{var}, '_',DT_runs{run}, '.eps');
            end
            disp(savename_eps)
            saveas(gcf,savename_eps,'epsc')
            %         export_fig(savename_eps,'-eps','-r600','-painters');
            close all
        end
    end
end
savename_struct = strcat(fig_path, '_', variable{var}, '_',DT_runs{1},'_',DT_runs{2}, '_plot_values.mat');
save(savename_struct, 'plot_values')
end
