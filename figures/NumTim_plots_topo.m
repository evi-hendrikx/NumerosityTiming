function NumTim_plots_topo(fig_path, NumTim_data, stat, DT_runs)
%% Make boxplot of the requested results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% fig_path: where you want to save the figure
% NumTim_data: made with the NumTim_load_data script. Contains all values
%       that you could need for all following analyses
% stat: structure containing statistical info from all the performed tests
% DT_runs: the two data runs you want to compare (options: TimingAll, TimingEven,
%       TimingOdd, NumerosityAll, NumerosityEven, NumerosityOdd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

n_fig = 0;
minAmountDataPoints = 1;

%% get information
[~,~,Run_1, Run_2] = NumTim_get_info(NumTim_data,DT_runs);
GNames = {{Run_1.Names{1:end-1}},{Run_2.Names{1:end-1}}};

if string(Run_1.Type) ~= string(Run_2.Type)
    new_runs = DT_runs;
else
    new_runs = {[Run_1.Type 'HalvesMean']};
end

for run = 1:length(new_runs)
   
    % prepare information
    parameters = []; data = [];

    if string(Run_1.Type) ~= string(Run_2.Type)
        data = {stat.data.(DT_runs{run}).topo.x0.x0};
    else
        ang1 = stat.data.(DT_runs{1}).topo.x0.x0;
        ang2 = stat.data.(DT_runs{2}).topo.x0.x0;

        % make sure means are not NaN if only 1 is NaN
        ang1(isnan(ang1)&~isnan(ang2)) = ang2(isnan(ang1)&~isnan(ang2));
        ang2(isnan(ang2)&~isnan(ang1)) = ang1(isnan(ang2)&~isnan(ang1));

        angles = rad2deg(circ_mean([deg2rad(ang1);deg2rad(ang2)]',[],2));
        if isnan(angles)
            angles = rad2deg(circ_mean([deg2rad(ang1),deg2rad(ang2)],[],2));
        end

        data = {angles};
    end
    parameters = [parameters {'x0x0'}];

    if isfield(stat.data.(DT_runs{run}).topo, 'y0') && string(Run_1.Type) ~= string(Run_2.Type)
        data = [data {stat.data.(DT_runs{run}).topo.y0.x0}];
        parameters = [parameters {'y0x0'}];
    end
    if isfield(stat.data.(DT_runs{run}).topo.x0, 'y0') && string(Run_1.Type) ~= string(Run_2.Type)
        data = [data {stat.data.(DT_runs{run}).topo.x0.y0}];
        parameters = [parameters {'x0y0'}];
    end
    if (isfield(stat.data.(DT_runs{run}).topo, 'y0') && isfield(stat.data.(DT_runs{run}).topo.x0, 'y0')) || isfield(stat.data.(DT_runs{run}).topo, 'y0') && isfield(stat.data.(DT_runs{run}).topo.y0, 'y0') 
        if string(Run_1.Type) ~= string(Run_2.Type)

            data = [data {stat.data.(DT_runs{run}).topo.y0.y0}];
        else
            ang1 = stat.data.(DT_runs{1}).topo.y0.y0;
            ang2 = stat.data.(DT_runs{2}).topo.y0.y0;

            % make sure means are not NaN if only 1 is NaN
            ang1(isnan(ang1)&~isnan(ang2)) = ang2(isnan(ang1)&~isnan(ang2));
            ang2(isnan(ang2)&~isnan(ang1)) = ang1(isnan(ang2)&~isnan(ang1));

            angles = rad2deg(circ_mean([deg2rad(ang1);deg2rad(ang2)]',[],2));
            if isnan(angles)
                angles = rad2deg(circ_mean([deg2rad(ang1),deg2rad(ang2)],[],2));
            end
            data = [data {angles}];

        end
        parameters = [parameters {'y0y0'}];
    end


    if string(Run_1.Type) == string(Run_2.Type)
        run = length(new_runs);
    end


    %% properties of figures
    circles_per_row = 3;
    size_figure = 14; % 90 single column, 140 1.5, 190 double
    size_subplots = (size_figure-2)/(circles_per_row*2);
    circles_per_column = 4;
    radius_circles = (size_subplots - 0.5)/2;

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

        %% determine properties subplots
        point_x=[];point_y=[];z_mean=[];point_y_median=[];

        % circular scatters
        colorMapScatter = linspace(1,10,8);

        left_data = data{par}(stat.data.(DT_runs{run}).hemi=='Left');
        left_maps = stat.data.(DT_runs{run}).map(stat.data.(DT_runs{run}).hemi=='Left');
        z_left = exp(1i*deg2rad(left_data)) * radius_circles;
        subj_left = stat.data.(DT_runs{run}).subj(stat.data.(DT_runs{run}).hemi=='Left');
        subNumsL = char(strip(subj_left,'S'));
        subjs_left = str2num(subNumsL(:));
        left_color = colorMapScatter(subjs_left);

        right_data = data{par}(stat.data.(DT_runs{run}).hemi=='Right');
        right_maps = stat.data.(DT_runs{run}).map(stat.data.(DT_runs{run}).hemi=='Right');
        z_right = exp(1i*deg2rad(right_data))* radius_circles;
        subj_right = stat.data.(DT_runs{run}).subj(stat.data.(DT_runs{run}).hemi=='Right');
        subNumsR = char(strip(subj_right,'S'));
        subjs_right = str2num(subNumsR(:));
        right_color = colorMapScatter(subjs_right);

        % summary stats with their CI
        nboot = 1000;
        for roi = 1:length(GNames{run})

            % same order as subjects later in color
            all_data_points = data{par}(stat.data.(DT_runs{run}).map == GNames{run}(roi));
            data_points_excl_NaN = all_data_points(~isnan(all_data_points));

            point_y = [point_y rad2deg(circ_mean(deg2rad(data_points_excl_NaN)))];
            if length(data_points_excl_NaN) > 1
                % "Interestingly, multiple quantities have been introduced as analogues to the linear standard deviation. First, the angular deviation is defined ass = sqrt(2 * (1 − R)). This quantity lies in the interval [0, sqrt(2)].
                % Alternatively, the circular standard deviation is defined as s0 = sqrt(−2 * ln(R)) and ranges from 0 to ∞.
                % Generally, the first measure is preferred, as it is bounded, but the two measures deviate little (Zar (1999))."
                one_stdev(roi) = rad2deg(circ_std(deg2rad(data_points_excl_NaN)));

                err_low = point_y(roi) - one_stdev(roi);
                if err_low < 0
                    err_low = 360 + err_low;
                end
                err_high = point_y(roi) + one_stdev(roi);
                if err_high > 360
                    err_high = err_high - 360;
                end

                error_bars(1,roi) = err_low; error_bars(2,roi) = err_high;

            else
                error_bars(:,roi) = [NaN NaN];
                one_stdev(roi) = NaN;
            end

            z_mean = [z_mean; exp(1i*deg2rad(point_y(roi))) * radius_circles];
            zz(:,roi) = exp(1i*linspace(0, 2*pi, 101))* radius_circles;


        end


        hold on

        circle_id = 0;

        for roi = 1:length(GNames{run})
            circle_id = circle_id + 1;

            sfh1 = subplot(circles_per_column,circles_per_row,circle_id);
            set(sfh1, 'Fontsize',7,'TickDir','out', 'FontName','Arial','Units','centimeters','Color','w','XColor','k','YColor','k','box','off');


            %% which subplot
            xposition_sp = 1 + mod(circle_id-1,3)*size_subplots;
            yposition_sp = 2+ circles_per_column * size_subplots - ceil(circle_id/3)*size_subplots;

            sfh1.Position = [xposition_sp yposition_sp size_subplots-0.1 size_subplots-0.1];

            cm = getPyPlot_cMap('Set2');
            colormap(gca, cm)
            set(sfh1,'CLim',[1,10]);

            hold on

            plot(real(zz(:,roi)), imag(zz(:,roi)), 'Color',[0.5, 0.5, 0.5])
            plot([-radius_circles-0.2 -radius_circles+0.2], [0 0], 'Color',[0.5, 0.5, 0.5],'LineStyle',':')
            plot([radius_circles-0.2 radius_circles+0.2], [0 0],  'Color',[0.5, 0.5, 0.5],'LineStyle',':')
            plot([0 0], [-radius_circles-0.2 -radius_circles+0.2], 'Color',[0.5, 0.5, 0.5],'LineStyle',':')
            plot([0 0], [radius_circles-0.2 radius_circles+0.2],  'Color',[0.5, 0.5, 0.5],'LineStyle',':');

            circd = @(radius,deg_ang)  [radius*cosd(deg_ang);  radius*sind(deg_ang)];       % Circle Function For Angles In Degrees
            N = 1000;                                                         % Number Of Points In Complete Circle
            if  one_stdev(roi) > 180 % 2 SDs spans the circle
                r_angl = linspace(error_bars(1,roi), error_bars(1,roi)+360, N);
            elseif error_bars(2,roi) <= error_bars(1,roi)
                r_angl = linspace(error_bars(1,roi), error_bars(2,roi)+360, N);                             % Angle Defining Arc Segment (radians)
            else
                r_angl = linspace(error_bars(1,roi), error_bars(2,roi), N);
            end
            xy_r = circd(radius_circles,r_angl);                                    % Matrix (2xN) Of (x,y) Coordinates
            p = plot(xy_r(1,:), xy_r(2,:), 'k-', 'LineWidth',1.5);
            p.Color = [0 0 0];

            scatter(real(z_left(left_maps== GNames{run}(roi))), imag(z_left(left_maps== GNames{run}(roi))), 10,left_color(left_maps== GNames{run}(roi)),'o','filled');
            scatter(real(z_right(right_maps== GNames{run}(roi))), imag(z_right(right_maps== GNames{run}(roi))), 10,right_color(right_maps== GNames{run}(roi)),'d','filled');
            scatter(real(z_mean(roi)), imag(z_mean(roi)), 30,'ks','filled'); %% make sure this is the right roi


            camroll(90)
            sfh1.XLabel.String = [GNames{run}(roi)];
            set(gca, 'XDir','reverse','xtick',[],'ytick',[])

            plot_values.(new_runs{run}).topo.(parameters{par}).(GNames{run}{roi}).errorbars = error_bars(:,roi);
            plot_values.(new_runs{run}).topo.(parameters{par}).(GNames{run}{roi}).r_angl = r_angl(:,roi);
            plot_values.(new_runs{run}).topo.(parameters{par}).(GNames{run}{roi}).mean = point_y(roi);
            plot_values.(new_runs{run}).topo.(parameters{par}).(GNames{run}{roi}).point = z_mean(roi);
            plot_values.(new_runs{run}).topo.(parameters{par}).(GNames{run}{roi}).SD = one_stdev(roi);



        end


        try
            savename_eps = strcat(fig_path, '_topo_',new_runs{run},'_',parameters{par}, '.eps');
        catch
            savename_eps = strcat(fig_path, '_topo_',new_runs{run}, '.eps');
        end
        disp(savename_eps)
        saveas(gcf,savename_eps,'epsc')
        %         export_fig(savename_eps,'-eps','-r600','-painters');
        close all
    end

    if string(Run_1.Type) == string(Run_2.Type)
        break
    end

end
savename_struct = strcat(fig_path, '_topo_',DT_runs{1},'_',DT_runs{2}, '_plot_values.mat');
save(savename_struct, 'plot_values')
end
