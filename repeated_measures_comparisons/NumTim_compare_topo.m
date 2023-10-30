function NumTim_compare_topo(fig_path, NumTim_data, whichCombi)
%% Compares correlations for Num-Time vs repeated measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% fig_path: path where you want figure to be saved
% NumTim_data: made with the NumTim_load_data script. Contains all values
%       that you could need for all following analyses
% whichCombi: 1 = compare to "All" maps; 2 = compare to most overlapping
%       map; 3 = compare to same map; 4 = compare to all overlapping maps;
%       5 = use voxel selection from other data structure (both use whichCombi = 4) 
%
% Output
% saves info in stat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minAmountDataPoints = 1;

%% To compute Bayes Factors, we need to compute average angles of NumEven-TimeOdd and NumOdd-TimeEven
% (this is already done for repeated measures wc=4 in create_table_for_R called
% by prepare_whichCombi

if whichCombi == 4
    stat_numEven_timeOdd = load('stat_NumerosityEven_TimingOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat');
    stat_numOdd_timeEven = load('stat_NumerosityOdd_TimingEven_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat');

    clear stat

    types = {'Numerosity', 'Timing'};
    runs = {'Even','Odd'};

    for type = 1:length(types)

        % gonna be the same for all because all compare numerosity (x0)
        % with timing (x0 and y0)
        par_1 = fieldnames(stat_numEven_timeOdd.stat.data.NumerosityEven.topo);
        par_2 = fieldnames(stat_numEven_timeOdd.stat.data.NumerosityEven.topo.(par_1{1}));

        for p1 = 1:length(par_1)

            for p2 = 1:length(par_2)

                rois = [];
                angs = [];
                hemi = [];
                subj = [];

                if type == 1
                    maps = unique(stat_numEven_timeOdd.stat.data.NumerosityEven.map);
                elseif type == 2
                    maps = unique(stat_numEven_timeOdd.stat.data.TimingOdd.map);
                end


                for map = 1:length(maps)
                    if type == 1
                        ang1 = stat_numEven_timeOdd.stat.data.NumerosityEven.topo.(par_1{p1}).(par_2{p2})(stat_numEven_timeOdd.stat.data.NumerosityEven.map == maps(map));
                        ang2 = stat_numOdd_timeEven.stat.data.NumerosityOdd.topo.(par_1{p1}).(par_2{p2})(stat_numOdd_timeEven.stat.data.NumerosityOdd.map == maps(map));
                    elseif type == 2
                        ang1 = stat_numEven_timeOdd.stat.data.TimingOdd.topo.(par_1{p1}).(par_2{p2})(stat_numEven_timeOdd.stat.data.TimingOdd.map == maps(map));
                        ang2 = stat_numOdd_timeEven.stat.data.TimingEven.topo.(par_1{p1}).(par_2{p2})(stat_numOdd_timeEven.stat.data.TimingEven.map == maps(map));
                    end

                    % make sure means are not NaN if only 1 is NaN
                    ang1(isnan(ang1)&~isnan(ang2)) = ang2(isnan(ang1)&~isnan(ang2));
                    ang2(isnan(ang2)&~isnan(ang1)) = ang1(isnan(ang2)&~isnan(ang1));

                    angles = rad2deg(circ_mean([deg2rad(ang1);deg2rad(ang2)]',[],2));

                    rois = [rois repmat(maps(map),1,length(angles))];
                    angs = [angs angles'];
                    if type == 1

                        hemi = [hemi stat_numEven_timeOdd.stat.data.NumerosityEven.hemi(stat_numEven_timeOdd.stat.data.NumerosityEven.map == maps(map))'];
                        subj = [subj stat_numEven_timeOdd.stat.data.NumerosityEven.subj(stat_numEven_timeOdd.stat.data.NumerosityEven.map == maps(map))'];
                    elseif type == 2
                        hemi = [hemi stat_numEven_timeOdd.stat.data.TimingOdd.hemi(stat_numEven_timeOdd.stat.data.TimingOdd.map == maps(map))'];
                        subj = [subj stat_numEven_timeOdd.stat.data.TimingOdd.subj(stat_numEven_timeOdd.stat.data.TimingOdd.map == maps(map))'];
                    end
                end
                DT_save = strcat(types{type},'Halves');

                stat.data.(DT_save).topo.(par_1{p1}).(par_2{p2}) = angs;
                stat.data.(DT_save).map = rois;

                stat.data.(DT_save).subj = subj;
                stat.data.(DT_save).hemi = hemi;

                tbl = table;
                tbl.map = rois';
                tbl.angles = angs';
                disp(tbl)
                if type == 1
                    par_tbl_path = strcat('stat_topoTbl', par_1{p1},'_',par_2{p2},'_meanRuns=1_NumerosityHalves_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.csv');
                elseif type == 2
                    par_tbl_path = strcat('stat_topoTbl', par_1{p1},'_',par_2{p2},'_meanRuns=1_TimingHalves_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.csv');
                end

                writetable(tbl,par_tbl_path)
            end
        end
    end
    savename_struct = strcat(strrep(fig_path,'fig','stat'),'.mat');
    save(savename_struct, 'stat')


elseif whichCombi == 5
    % Num odd and num even in data loaded here give similar topo values, because done in
    % the same voxel selection
    num_stat_numEven_timeOdd = load('stat_NumerosityEven_NumerosityOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=5_selection=NumerosityEven_TimingOdd_topo_measure=mean.mat');
    num_stat_numOdd_timeEven = load('stat_NumerosityEven_NumerosityOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=5_selection=NumerosityOdd_TimingEven_topo_measure=mean.mat');
    % time odd and time even here give similar topo values, because done in
    % the same voxel selection
    time_stat_numEven_timeOdd = load('stat_TimingEven_TimingOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=5_selection=NumerosityEven_TimingOdd_topo_measure=mean.mat');
    time_stat_numOdd_timeEven = load('stat_TimingEven_TimingOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=5_selection=NumerosityOdd_TimingEven_topo_measure=mean.mat');


    clear stat
    
    % average repeated measures over different selections (data type here
    % arbitrary, as long as it's in the repeated measures, see comment
    % above)
    num_NETO_x0x0 = num_stat_numEven_timeOdd.stat.data.NumerosityEven.topo.x0.x0;
    num_NOTE_x0x0 = num_stat_numOdd_timeEven.stat.data.NumerosityOdd.topo.x0.x0;
    num_NETO_x0x0(isnan(num_NETO_x0x0)& ~isnan(num_NOTE_x0x0)) = num_NOTE_x0x0(isnan(num_NETO_x0x0)& ~isnan(num_NOTE_x0x0));
    num_NOTE_x0x0(~isnan(num_NETO_x0x0)& isnan(num_NOTE_x0x0)) = num_NETO_x0x0(~isnan(num_NETO_x0x0)& isnan(num_NOTE_x0x0));
    rm_num_topo_x0x0 = rad2deg(circ_mean([deg2rad(num_NETO_x0x0);deg2rad(num_NOTE_x0x0)]',[],2));

    tbl = table;
    tbl.map = num_stat_numEven_timeOdd.stat.data.NumerosityEven.map; % same order for all
    tbl.angles = rm_num_topo_x0x0;
    disp(tbl)
    par_tbl_path = strcat('stat_topoTblx0x0_meanRuns=1_NumerosityEven_NumerosityOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=5_topo_measure=mean.csv');
    writetable(tbl,par_tbl_path)
    stat.data.NumerosityRM.topo.x0.x0 = rm_num_topo_x0x0;

    time_NETO_x0x0 = time_stat_numEven_timeOdd.stat.data.TimingOdd.topo.x0.x0;
    time_NOTE_x0x0 = time_stat_numOdd_timeEven.stat.data.TimingEven.topo.x0.x0;
    time_NETO_x0x0(isnan(time_NETO_x0x0)& ~isnan(time_NOTE_x0x0)) = time_NOTE_x0x0(isnan(time_NETO_x0x0)& ~isnan(time_NOTE_x0x0));
    time_NOTE_x0x0(~isnan(time_NETO_x0x0)& isnan(time_NOTE_x0x0)) = time_NETO_x0x0(~isnan(time_NETO_x0x0)& isnan(time_NOTE_x0x0));
    rm_time_topo_x0x0 = rad2deg(circ_mean([deg2rad(time_NETO_x0x0);deg2rad(time_NOTE_x0x0)]',[],2));

    tbl = table;
    tbl.map = time_stat_numEven_timeOdd.stat.data.TimingOdd.map; % same order for all
    tbl.angles = rm_time_topo_x0x0;
    disp(tbl)
    par_tbl_path = strcat('stat_topoTblx0x0_meanRuns=1_TimingEven_TimingOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=5_topo_measure=mean.csv');
    writetable(tbl,par_tbl_path)
    stat.data.TimingRM.topo.x0.x0 = rm_time_topo_x0x0;


    time_NETO_y0y0 = time_stat_numEven_timeOdd.stat.data.TimingOdd.topo.y0.y0;
    time_NOTE_y0y0 = time_stat_numOdd_timeEven.stat.data.TimingEven.topo.y0.y0;
    time_NETO_y0y0(isnan(time_NETO_y0y0)& ~isnan(time_NOTE_y0y0)) = time_NOTE_y0y0(isnan(time_NETO_y0y0)& ~isnan(time_NOTE_y0y0));
    time_NOTE_y0y0(~isnan(time_NETO_y0y0)& isnan(time_NOTE_y0y0)) = time_NETO_y0y0(~isnan(time_NETO_y0y0)& isnan(time_NOTE_y0y0));
    rm_time_topo_y0y0 = rad2deg(circ_mean([deg2rad(time_NETO_y0y0);deg2rad(time_NOTE_y0y0)]',[],2));

    tbl = table;
    tbl.map = time_stat_numEven_timeOdd.stat.data.TimingOdd.map; % same order for all
    tbl.angles = rm_time_topo_y0y0;
    disp(tbl)
    par_tbl_path = strcat('stat_topoTbly0y0_meanRuns=1_TimingEven_TimingOdd_minVEtime=0.2_minVEnum=0.3_whichCombi=5_topo_measure=mean.csv');
    writetable(tbl,par_tbl_path)
    stat.data.TimingRM.topo.y0.y0 = rm_time_topo_y0y0;


    % make averages different runs (selection criteria that were used)
    num_maps = num_stat_numOdd_timeEven.stat.data.NumerosityEven.map;
    time_maps = time_stat_numOdd_timeEven.stat.data.TimingEven.map;
    
    num_NETO_selection_x0x0 = num_stat_numEven_timeOdd.stat.data.NumerosityEven.selection_topo.x0.x0;
    num_NOTE_selection_x0x0 = num_stat_numOdd_timeEven.stat.data.NumerosityOdd.selection_topo.x0.x0;
    num_NETO_selection_x0x0(isnan(num_NETO_selection_x0x0)& ~isnan(num_NOTE_selection_x0x0)) = num_NOTE_selection_x0x0(isnan(num_NETO_selection_x0x0)& ~isnan(num_NOTE_selection_x0x0));
    num_NOTE_selection_x0x0(~isnan(num_NETO_selection_x0x0)& isnan(num_NOTE_selection_x0x0)) = num_NETO_selection_x0x0(~isnan(num_NETO_selection_x0x0)& isnan(num_NOTE_selection_x0x0));
    diff_topo_x0x0_num = rad2deg(circ_mean([deg2rad(num_NETO_selection_x0x0);deg2rad(num_NOTE_selection_x0x0)]',[],2));

    tbl = table;
    tbl.map = num_maps;
    tbl.angles = diff_topo_x0x0_num;
    disp(tbl)
    par_tbl_path = strcat('stat_topoTblx0x0_meanRuns=1_NumerosityHalves_minVEtime=0.2_minVEnum=0.3_whichCombi=5_topo_measure=mean.csv');
    writetable(tbl,par_tbl_path)
    stat.data.NumerosityHalves.topo.x0.x0 = diff_topo_x0x0_num;

    
    num_NETO_selection_x0y0 = num_stat_numEven_timeOdd.stat.data.NumerosityEven.selection_topo.x0.y0;
    num_NOTE_selection_x0y0 = num_stat_numOdd_timeEven.stat.data.NumerosityOdd.selection_topo.x0.y0;
    num_NETO_selection_x0y0(isnan(num_NETO_selection_x0y0)& ~isnan(num_NOTE_selection_x0y0)) = num_NOTE_selection_x0y0(isnan(num_NETO_selection_x0y0)& ~isnan(num_NOTE_selection_x0y0));
    num_NOTE_selection_x0y0(~isnan(num_NETO_selection_x0y0)& isnan(num_NOTE_selection_x0y0)) = num_NETO_selection_x0y0(~isnan(num_NETO_selection_x0y0)& isnan(num_NOTE_selection_x0y0));
    diff_topo_x0y0_num = rad2deg(circ_mean([deg2rad(num_NETO_selection_x0y0);deg2rad(num_NOTE_selection_x0y0)]',[],2));

    tbl = table;
    tbl.map = num_maps;
    tbl.angles = diff_topo_x0y0_num;
    disp(tbl)
    par_tbl_path = strcat('stat_topoTblx0y0_meanRuns=1_NumerosityHalves_minVEtime=0.2_minVEnum=0.3_whichCombi=5_topo_measure=mean.csv');
    writetable(tbl,par_tbl_path)
    stat.data.NumerosityHalves.topo.x0.y0 = diff_topo_x0y0_num;

    
    time_NETO_selection_x0x0 = time_stat_numEven_timeOdd.stat.data.TimingOdd.selection_topo.x0.x0;
    time_NOTE_selection_x0x0 = time_stat_numOdd_timeEven.stat.data.TimingEven.selection_topo.x0.x0;
    time_NETO_selection_x0x0(isnan(time_NETO_selection_x0x0)& ~isnan(time_NOTE_selection_x0x0)) = time_NOTE_selection_x0x0(isnan(time_NETO_selection_x0x0)& ~isnan(time_NOTE_selection_x0x0));
    time_NOTE_selection_x0x0(~isnan(time_NETO_selection_x0x0)& isnan(time_NOTE_selection_x0x0)) = time_NETO_selection_x0x0(~isnan(time_NETO_selection_x0x0)& isnan(time_NOTE_selection_x0x0));
    diff_topo_x0x0_time = rad2deg(circ_mean([deg2rad(time_NETO_selection_x0x0);deg2rad(time_NOTE_selection_x0x0)]',[],2));

    tbl = table;
    tbl.map = time_maps;
    tbl.angles = diff_topo_x0x0_time;
    disp(tbl)
    par_tbl_path = strcat('stat_topoTblx0x0_meanRuns=1_TimingHalves_minVEtime=0.2_minVEnum=0.3_whichCombi=5_topo_measure=mean.csv');
    writetable(tbl,par_tbl_path)
    stat.data.TimingHalves.topo.x0.x0 = diff_topo_x0x0_time;
    
    time_NETO_selection_x0y0 = time_stat_numEven_timeOdd.stat.data.TimingOdd.selection_topo.x0.y0;
    time_NOTE_selection_x0y0 = time_stat_numOdd_timeEven.stat.data.TimingEven.selection_topo.x0.y0;
    time_NETO_selection_x0y0(isnan(time_NETO_selection_x0y0)& ~isnan(time_NOTE_selection_x0y0)) = time_NOTE_selection_x0y0(isnan(time_NETO_selection_x0y0)& ~isnan(time_NOTE_selection_x0y0));
    time_NOTE_selection_x0y0(~isnan(time_NETO_selection_x0y0)& isnan(time_NOTE_selection_x0y0)) = time_NETO_selection_x0y0(~isnan(time_NETO_selection_x0y0)& isnan(time_NOTE_selection_x0y0));
    diff_topo_x0y0_time = rad2deg(circ_mean([deg2rad(time_NETO_selection_x0y0);deg2rad(time_NOTE_selection_x0y0)]',[],2));

    tbl = table;
    tbl.map = time_maps;
    tbl.angles = diff_topo_x0y0_time;
    disp(tbl)
    par_tbl_path = strcat('stat_topoTblx0y0_meanRuns=1_TimingHalves_minVEtime=0.2_minVEnum=0.3_whichCombi=5_topo_measure=mean.csv');
    writetable(tbl,par_tbl_path)
    stat.data.TimingHalves.topo.x0.y0 = diff_topo_x0y0_time;
    

    % build stat for plot
 
    stat.data.TimingHalves.map = time_maps;
    stat.data.NumerosityHalves.map = num_maps;
    stat.data.TimingHalves.subj = time_stat_numEven_timeOdd.stat.data.TimingOdd.subj;
    stat.data.NumerosityHalves.subj = num_stat_numEven_timeOdd.stat.data.NumerosityEven.subj;
    stat.data.TimingHalves.hemi = time_stat_numEven_timeOdd.stat.data.TimingOdd.hemi;
    stat.data.NumerosityHalves.hemi = num_stat_numEven_timeOdd.stat.data.NumerosityEven.hemi;

    stat.data.TimingRM.map = time_maps;
    stat.data.NumerosityRM.map = num_maps;
    stat.data.TimingRM.subj = time_stat_numEven_timeOdd.stat.data.TimingOdd.subj;
    stat.data.NumerosityRM.subj = num_stat_numEven_timeOdd.stat.data.NumerosityEven.subj;
    stat.data.TimingRM.hemi = time_stat_numEven_timeOdd.stat.data.TimingOdd.hemi;
    stat.data.NumerosityRM.hemi = num_stat_numEven_timeOdd.stat.data.NumerosityEven.hemi;
    
    savename_struct = strcat(strrep(fig_path,'fig','stat'),'.mat');
    save(savename_struct, 'stat')

else
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN STATISTICS IN R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BFs different from previous (hard coded (sorry)) BFs in other scripts,
% because these here concern BFs over halves of the data (and previously
% everything together)

BF_info = readmatrix(strcat('BF_wc',num2str(whichCombi),'_x0x0_numerosityHalves.csv'), NumHeaderLines = 1, OutputType = 'string');

BF_roi_names = BF_info(:,2);
BF_factors = BF_info(:,end);
BF_factors = cellfun(@str2num,BF_factors);

BF_rois.NumerosityHalves = BF_roi_names;
BF.NumerosityHalves.x0x0 = BF_factors;


BF_info = readmatrix(strcat('BF_wc',num2str(whichCombi),'_x0y0_numerosityHalves.csv'), NumHeaderLines = 1, OutputType = 'string');
BF_factors = BF_info(:,end);
BF_factors = cellfun(@str2num,BF_factors);

BF.NumerosityHalves.x0y0 = BF_factors;

BF_info = readmatrix(strcat('BF_wc',num2str(whichCombi),'_x0x0_timingHalves.csv'), NumHeaderLines = 1, OutputType = 'string');
BF_roi_names = BF_info(:,2);
BF_factors = BF_info(:,end);
BF_factors = cellfun(@str2num,BF_factors);

BF_rois.TimingHalves = BF_roi_names;
BF.TimingHalves.x0x0 = BF_factors;

BF_info = readmatrix(strcat('BF_wc',num2str(whichCombi),'_x0y0_timingHalves.csv'), NumHeaderLines = 1, OutputType = 'string');
BF_factors = BF_info(:,end);
BF_factors = cellfun(@str2num,BF_factors);

BF.TimingHalves.x0y0 = BF_factors;

% repeated measures
BF_info = readmatrix(strcat('BF_wc',num2str(whichCombi),'_x0x0_rm_numerosity.csv'), NumHeaderLines = 1, OutputType = 'string');
BF_roi_names = BF_info(:,2);
BF_factors = BF_info(:,end);
BF_factors = cellfun(@str2num,BF_factors);

[~,~,index_rois]=intersect(BF_rois.NumerosityHalves,BF_roi_names,'stable');

BF_rm.NumerosityHalves.x0x0 = BF_factors(index_rois);

BF_info = readmatrix(strcat('BF_wc',num2str(whichCombi),'_x0x0_rm_timing.csv'), NumHeaderLines = 1, OutputType = 'string');
BF_roi_names = BF_info(:,2);
BF_factors = BF_info(:,end);
BF_factors = cellfun(@str2num,BF_factors);

clear index_rois
[~,~,index_rois]=intersect(BF_rois.TimingHalves,BF_roi_names,'stable');


BF_rm.TimingHalves.x0x0 = BF_factors(index_rois);

BF_info = readmatrix(strcat('BF_wc',num2str(whichCombi),'_y0y0_rm_timing.csv'), NumHeaderLines = 1, OutputType = 'string');
BF_roi_names = BF_info(:,2);
BF_factors = BF_info(:,end);
BF_factors = cellfun(@str2num,BF_factors);

clear index_rois
[~,~,index_rois]=intersect(BF_rois.TimingHalves,BF_roi_names,'stable');

BF_rm.TimingHalves.y0y0 = BF_factors(index_rois);





%% check normality
% Num-Dur and Num-Per
DT_runs = {'NumerosityHalves', 'TimingHalves'};
diff_parameters = {'x0x0','x0y0'};
NumerosityHalves_parameters = {'x0x0', 'x0x0'};
TimingHalves_parameters = {'x0x0','y0y0'};
p_norms = [];

for run = 1:length(DT_runs)
    for par = 1:length(diff_parameters)
        eval(['same_parameters =', DT_runs{run}, '_parameters;'])

        diff_BFs = BF.(DT_runs{run}).(diff_parameters{par});
        same_BFs = BF_rm.(DT_runs{run}).(same_parameters{par});

        try
            [~,stat.(DT_runs{run}).topo.(diff_parameters{par}).(same_parameters{par}).normality,~] = force_swtest(same_BFs - diff_BFs,0.05,1);
        catch
            stat.(DT_runs{run}).topo.(diff_parameters{par}).(same_parameters{par}).normality = 0;
        end
        p_norms = [p_norms stat.(DT_runs{run}).topo.(diff_parameters{par}).(same_parameters{par}).normality];


    end
end
disp(p_norms)

%% comparisons
nboot = 1000;

for run = 1:length(DT_runs)

    % FDR correction takes place within numerosity maps when
    % NumEven-DurEven vs. NumEven-NumOdd and NumEven-PerEven vs. NumEven-NumOdd
    % NOT INCLUDED THERE:
    % NumOdd-DurOdd and NumOdd-PerOdd, since these are done on
    % different voxel selections (overlap with NumEven-TimeEven and
    % NumOdd-TimeOdd, respectively)
    pvals= [];
    store_rois = [];
    comparison = [];
    comp = [];

    for par = 1:length(diff_parameters)
        eval(['same_parameters =', DT_runs{run}, '_parameters;'])

        diff_BFs = BF.(DT_runs{run}).(diff_parameters{par});
        same_BFs = BF_rm.(DT_runs{run}).(same_parameters{par});

        stat.(DT_runs{run}).topo.(diff_parameters{par}).(same_parameters{par}).data_diff = diff_BFs;
        stat.(DT_runs{run}).topo.(diff_parameters{par}).(same_parameters{par}).data_rm = same_BFs;

        stat.(DT_runs{run}).topo.(diff_parameters{par}).(same_parameters{par}).median_diff = median(diff_BFs);
        stat.(DT_runs{run}).topo.(diff_parameters{par}).(same_parameters{par}).median_rm = median(same_BFs);
        stat.(DT_runs{run}).topo.(diff_parameters{par}).(same_parameters{par}).CI_diff = bootci(nboot, {@nanmedian, diff_BFs},'alpha',0.05);
        stat.(DT_runs{run}).topo.(diff_parameters{par}).(same_parameters{par}).CI_rm = bootci(nboot, {@nanmedian, same_BFs},'alpha',0.05);

        [p,stat.(DT_runs{run}).topo.(diff_parameters{par}).(same_parameters{par}).h,stat.(DT_runs{run}).topo.(diff_parameters{par}).(same_parameters{par}).stats] = signrank(same_BFs,diff_BFs,'method','approximate');

        stat.(DT_runs{run}).topo.(diff_parameters{par}).(same_parameters{par}).p = p;

        pvals = [pvals p];

    end

    adj_p = NaN(1,length(pvals));
    [~,~,~,adj_p(~isnan(pvals))] = fdr_bh(pvals(~isnan(pvals)));

    for par = 1:length(diff_parameters)
        eval(['same_parameters =', DT_runs{run}, '_parameters;'])
        stat.(DT_runs{run}).topo.(diff_parameters{par}).(same_parameters{par}).adj_p = adj_p(par);
    end

end

save(savename_struct, 'stat')

%% plots for angels over halves
NumTim_plots_topo(fig_path, NumTim_data, stat, {'NumerosityHalves','TimingHalves'})
if whichCombi == 5
    NumTim_plots_topo(fig_path, NumTim_data, stat, {'NumerosityRM','NumerosityRM'})
    NumTim_plots_topo(fig_path, NumTim_data, stat, {'TimingRM','TimingRM'})
end

