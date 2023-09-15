function stat = NumTim_chance_anova_stats(stat_path, whichCombi,stat,DT_runs)
%% Calculates stats for provided information (proportion, correlation & topographic angles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% stat_path: path where you want stat to be saved
% whichCombi: 1 = compare to "All" maps; 2 = compare to most overlapping
%       map; 3 = compare to same map; 4 = compare to all overlapping maps
% stat: structure containing statistical info currently containing data on
%       which tests should be applies
% DT_runs: the two data runs you want to compare (options: TimingAll, TimingEven,
%       TimingOdd, NumerosityAll, NumerosityEven, NumerosityOdd)
%
% Output
% stat: structure containing statistical info from all the performed tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check if it exists
if (exist(stat_path, 'file') == 2)
    load(stat_path)
    return
end

% general info
minAmountDataPoints = 1;
variable = {'proportion','correlation'};


for var = 1:length(variable)

    %% conveniently store data
    for run = 1:length(DT_runs)
        eval(strcat('parameters',char(num2str(run)),' = []; data',char(num2str(run)),' = [];'))
        if variable{var} == "proportion"
            if whichCombi == 1 || whichCombi == 4
                eval(strcat('chance_data',char(num2str(run)),' = {stat.data.(DT_runs{run}).chance_proportion};'))
            end
            eval(strcat('data',char(num2str(run)),' = [{stat.data.(DT_runs{run}).prop}];'))
            eval(strcat('size_data',char(num2str(run)),' = {stat.data.(DT_runs{run}).size};'))
            eval(strcat('parameters',char(num2str(run)),' = [{char("proportion")}];'))

            [roi_names.(DT_runs{run}), ~, roi_id_of_map.(DT_runs{run})] = unique(stat.data.(DT_runs{run}).map,'stable');

        elseif variable{var} == "correlation"
            eval(strcat('data',char(num2str(run)),'= {stat.data.(DT_runs{run}).correlation.x0.x0};'))
            eval(strcat('parameters',char(num2str(run)),' = [parameters',char(num2str(run)),' {char("x0x0")}];'))
            if isfield(stat.data.(DT_runs{run}).correlation, 'y0') && ~(contains(DT_runs{run},"Timing") && contains(DT_runs{3-run},"Timing"))
                eval(strcat('data',char(num2str(run)),' = [data',char(num2str(run)),' {stat.data.(DT_runs{run}).correlation.y0.x0}];'))
                eval(strcat('parameters',char(num2str(run)),' = [parameters',char(num2str(run)),' {char("y0x0")}];'))
            end
            if isfield(stat.data.(DT_runs{run}).correlation.x0, 'y0') && ~(contains(DT_runs{run},"Timing") && contains(DT_runs{3-run},"Timing"))
                eval(strcat('data',char(num2str(run)),' = [data',char(num2str(run)),' {stat.data.(DT_runs{run}).correlation.x0.y0}];'))
                eval(strcat('parameters',char(num2str(run)),' = [parameters',char(num2str(run)),' {char("x0y0")}];'))
            end
            if isfield(stat.data.(DT_runs{run}).correlation, 'y0') && isfield(stat.data.(DT_runs{run}).correlation.x0, 'y0')
                eval(strcat('data',char(num2str(run)),' = [data',char(num2str(run)),' {stat.data.(DT_runs{run}).correlation.y0.y0}];'))
                eval(strcat('parameters',char(num2str(run)),' = [parameters',char(num2str(run)),' {char("y0y0")}];'))
            end

            % kick out areas with only 1 datapoint from correlational
            % analyses
            [roi_names.(DT_runs{run}), ~, roi_id_of_map.(DT_runs{run})] = unique(stat.data.(DT_runs{run}).map,'stable');
            loop_roi_names = roi_names.(DT_runs{run}); loop_roi_ids = roi_id_of_map.(DT_runs{run});
            kicked_rois = 0;
            for roi = 1:length(loop_roi_names)
                for par = 1:length(eval(['parameters',char(num2str(run))]))
                    if eval(['sum(~isnan(data',char(num2str(run)),'{par}(loop_roi_ids == roi))) <= minAmountDataPoints'])
                        eval(['data',char(num2str(run)),'{par}(loop_roi_ids == roi) = NaN;'])
                        % only do this once, if missing for x0x0, will also
                        % miss for other correlations and vice versa
                        if par == 1
                            roi_names.(DT_runs{run})(roi - kicked_rois) = [];
                            roi_id_of_map.(DT_runs{run})(loop_roi_ids == roi) = NaN;
                            roi_id_of_map.(DT_runs{run})(loop_roi_ids > roi) = roi_id_of_map.(DT_runs{run})(loop_roi_ids > roi) - 1;
                            kicked_rois = kicked_rois + 1;
                        end
                    end
                end
            end
        end
    end

    %% normality tests chance distribution
    pvals_norm_chance = [];
    for run = 1:length(DT_runs)

        % what to test against
        % chance distribution doesn't conceptually make sense if
        % you're not looking at "All" overlap for proportions
        if variable{var} == "proportion" && whichCombi ~= 1 && whichCombi ~= 4
            continue
        end

        current_roi_names = roi_names.(DT_runs{run});
        current_roi_id_of_map = roi_id_of_map.(DT_runs{run});

        parameters = eval(strcat('parameters',char(num2str(run))));
        data = eval(strcat('data',char(num2str(run))));

        for par = 1:eval(strcat('length(parameters',char(num2str(run)),')'))
            % rois info
            if variable{var} == "proportion"
                chance_distribution = strcat('chance_data',char(num2str(run)),'{par}(current_roi_id_of_map == roi)');
                chance_distribution_all = strcat('chance_data',char(num2str(run)),'{par}');
            else % correlations
                chance_distribution = '0';
                chance_distribution_all = '0';
            end

            % check normality for all maps together
            eval(strcat('[~,stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).All.p_norm,~] = force_swtest(data{par}-',chance_distribution_all,',0.05,1);'))

            % check normality per map
            for roi = 1:length(current_roi_names)
                try
                    eval(strcat('[~,stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).(current_roi_names{roi}).p_norm,~] = force_swtest(data{par}(current_roi_id_of_map == roi)-',chance_distribution,',0.05,1);'))
                catch
                    stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).(current_roi_names{roi}).p_norm = 0;
                end

                % store them just to easily check whether <0.05 or not
                pvals_norm_chance = [pvals_norm_chance stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).(current_roi_names{roi}).p_norm];
            end

        end
    end

    %% statistics against chance
    for run = 1:length(DT_runs)
        % what to test against
        % chance distribution doesn't conceptually make sense if
        % you're not looking at "All" overlap for proportions
        if variable{var} == "proportion" && whichCombi ~= 1 && whichCombi ~= 4
            continue
        end

        parameters = eval(strcat('parameters',char(num2str(run))));
        data = eval(strcat('data',char(num2str(run))));
        current_roi_names = roi_names.(DT_runs{run});
        current_roi_id_of_map = roi_id_of_map.(DT_runs{run});

        % done for duration AND period in same maps in order to apply FDR
        pvals_chance = [];
        pvals_chance_all = [];

        for par = 1:length(parameters)

            % tests against 0 correlation
            % or check whether proportion is significantly better than chance
            % chance distribution doesn't conceptually make sense if
            % you're not looking at "All" overlap for proportions
            if parameters{par} == "proportion" % proportion
                chance_distribution = strcat('chance_data',char(num2str(run)),'{par}(current_roi_id_of_map == roi)');
                chance_distribution_all = strcat('chance_data',char(num2str(run)),'{par}');
                tail = char("'right'");
            else
                chance_distribution = '0';
                chance_distribution_all = '0';
                tail = char("'both'");
            end


            % normally distributed
            if isempty(find(pvals_norm_chance(~isnan(pvals_norm_chance))<0.05, 1))
                results_stats = '[~,stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).(current_roi_names{roi}).p,stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).(current_roi_names{roi}).ci,stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).(current_roi_names{roi}).stats]';
                results_stats_all = '[~,stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).All.p,stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).All.ci,stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).All.stats]';
                test ='ttest';
                method = [];

                % not normally distributed
            else
                test ='signrank';
                results_stats = '[stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).(current_roi_names{roi}).p,stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).(current_roi_names{roi}).h,stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).(current_roi_names{roi}).stats]';
                results_stats_all = '[stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).All.p,stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).All.h,stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).All.stats]';
                method = [char(",'method',"),char("'approximate'")];
            end

            eval(strcat(results_stats_all,'=',test,'(data{par},',chance_distribution_all,char(",'tail',"),tail,method,');'))
            pvals_chance_all = [pvals_chance_all stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).All.p];

            for roi = 1:length(current_roi_names)
                eval(strcat(results_stats,'=',test,'(data{par}(current_roi_id_of_map == roi),',chance_distribution,char(",'tail',"),tail,method,');'))
                pvals_chance = [pvals_chance stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).(current_roi_names{roi}).p];
            end

        end

        adj_p_all = nan(size(pvals_chance_all));
        [~,~,~,adj_p_all(~isnan(pvals_chance_all))] = fdr_bh(pvals_chance_all(~isnan(pvals_chance_all)));

        adj_p = nan(size(pvals_chance));
        [~,~,~,adj_p(~isnan(pvals_chance))] = fdr_bh(pvals_chance(~isnan(pvals_chance)));

        % done for duration AND period in same maps in order to apply FDR
        for par = 1:length(parameters)
            stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).All.adj_p = adj_p_all(par);
            for roi = 1:length(current_roi_names)
                stat.different_from_chance.(DT_runs{run}).(variable{var}).(parameters{par}).(current_roi_names{roi}).adj_p = adj_p(roi + length(current_roi_names)*(par-1));
            end
        end
    end


    p_norms_ANOVA = [];
    terms(1,1) = 1; % main effect subject
    terms(2,2) = 1; % main effect map
    terms(3,3) = 1; % main effect hemisphere
    %% statistics between maps: normality
    for run = 1:length(DT_runs)

        parameters = eval(strcat('parameters',char(num2str(run))));
        data = eval(strcat('data',char(num2str(run))));

        for par = 1:length(parameters)
            if variable{var} == "proportion" && parameters{par} == "proportion"
                store = '(variable{var})';
            else
                store = '(variable{var}).(parameters{par})';
            end

            %% ANOVA including post hocs: are the hemispheres/ maps different from each other
            % normality test on the residuals
            [~, tbl, stats, ~] = anovan(data{par}, {stat.data.(DT_runs{run}).map, stat.data.(DT_runs{run}).hemi, stat.data.(DT_runs{run}).subj}, 'varnames',{'map','hemisphere','subject'});
            [~,p_norm,~] = force_swtest(stats.resid,0.05,1);

            p_norms_ANOVA = [p_norms_ANOVA p_norm];


            eval(strcat('stat.ANOVA.(DT_runs{run}).', store,'.tbl = tbl;'))
            eval(strcat('stat.ANOVA.(DT_runs{run}).', store,'.stats = stats;'))
            eval(strcat('stat.ANOVA.(DT_runs{run}).', store,'.terms = terms;'))
            eval(strcat('stat.ANOVA.(DT_runs{run}).', store,'.p_norm = p_norm;'))
        end
    end

    for run = 1:length(DT_runs)
        parameters = eval(strcat('parameters',char(num2str(run))));
        data = eval(strcat('data',char(num2str(run))));
        current_roi_names = roi_names.(DT_runs{run});
        current_roi_id_of_map = roi_id_of_map.(DT_runs{run});

        for par = 1:length(parameters)

            if variable{var} == "proportion" && parameters{par} == "proportion"
                store = '(variable{var})';
            else
                store = '(variable{var}).(parameters{par})';
            end

            % post hocs
            if eval(strcat('stat.ANOVA.(DT_runs{run}).',store,'.tbl{string(stat.ANOVA.(DT_runs{run}).',store,'.tbl(:,1))=="map",string(stat.ANOVA.(DT_runs{run}).',store,'.tbl(1,:))=="Prob>F"} < 0.05'))
                if isempty(find(p_norms_ANOVA<0.05, 1))
                    [results,~,~,gnames] = multcompare(eval(strcat('stat.ANOVA.(DT_runs{run}).', store,'.stats')),'Dimension',[1]);
                    eval(strcat('stat.posthoc.(DT_runs{run}).',store,'.gnames = gnames;'))
                    eval(strcat('stat.posthoc.(DT_runs{run}).',store,'.multcomp = results;'))
                    close all
                else
                    dunn_roi_id_of_map  = current_roi_id_of_map(~isnan(data{par}));
                    dunn_data = data{par}(~isnan(data{par}));

                    [roi_id_of_map_grouped, grouping_id_maps] = sort(dunn_roi_id_of_map);
                    dunn_data = reshape(dunn_data(grouping_id_maps),1,length(dunn_data));

                    %roi_id_of_map_grouped = roi_id_of_map_grouped(~isnan(dunn_data));
                    roi_id_of_map_grouped = reshape(roi_id_of_map_grouped,1,length(dunn_data));


                    dunns = dunn(dunn_data,roi_id_of_map_grouped);

                    eval(strcat('stat.posthoc.(DT_runs{run}).',store,'.roi_names = current_roi_names;'))
                    eval(strcat('stat.posthoc.(DT_runs{run}).', store,'.tbl = dunns;'))
                end
            end
        end
    end

    %% calculate relation proportion and size
    if variable{var} == "proportion"
        pvals_subj = []; pvals_hemi = []; pvals_size = [];
        for run = 1:length(DT_runs)

            not_nan = ~isnan(stat.data.(DT_runs{run}).prop);

            covariates = {stat.data.(DT_runs{run}).size, stat.data.(DT_runs{run}).hemi,stat.data.(DT_runs{run}).subj};
            covariate_names = {'size', 'hemi', 'subj'};

            [~,~,stat.different_from_chance.(DT_runs{run}).(variable{var}).size_relation.hemi_ctab,~] = aoctool(stat.data.(DT_runs{run}).size(not_nan),stat.data.(DT_runs{run}).prop(not_nan),stat.data.(DT_runs{run}).hemi(not_nan), 0.05, 'size','prop','hemi','on', 'parallel lines');
            [~,~,stat.different_from_chance.(DT_runs{run}).(variable{var}).size_relation.subj_ctab,~] = aoctool(stat.data.(DT_runs{run}).size(not_nan),stat.data.(DT_runs{run}).prop(not_nan),stat.data.(DT_runs{run}).subj(not_nan), 0.05, 'size','prop','subj','on', 'parallel lines');

            [~, tbl, stats, ~] = anovan(stat.data.(DT_runs{run}).prop, covariates,'Continuous',1,'varnames',covariate_names);
            stat.different_from_chance.(DT_runs{run}).(variable{var}).size_relation.tbl = tbl;
            stat.different_from_chance.(DT_runs{run}).(variable{var}).size_relation.stats = stats;

            pvals_subj = [pvals_subj stat.different_from_chance.(DT_runs{run}).(variable{var}).size_relation.tbl{(string(stat.different_from_chance.(DT_runs{run}).(variable{var}).size_relation.tbl(:,1))=="subj"),7}];
            pvals_hemi = [pvals_hemi stat.different_from_chance.(DT_runs{run}).(variable{var}).size_relation.tbl{(string(stat.different_from_chance.(DT_runs{run}).(variable{var}).size_relation.tbl(:,1))=="hemi"),7}];
            pvals_size = [pvals_size stat.different_from_chance.(DT_runs{run}).(variable{var}).size_relation.tbl{(string(stat.different_from_chance.(DT_runs{run}).(variable{var}).size_relation.tbl(:,1))=="size"),7}];

            close all
        end

        [~,~,~,adjp_subj] = fdr_bh(pvals_subj);
        [~,~,~,adjp_hemi] = fdr_bh(pvals_hemi);
        [~,~,~,adjp_size] = fdr_bh(pvals_size);

        for run = 1:length(DT_runs)

        stat.different_from_chance.(DT_runs{run}).(variable{var}).size_relation.adj_p.subj = adjp_subj(run);
        stat.different_from_chance.(DT_runs{run}).(variable{var}).size_relation.adj_p.hemi = adjp_hemi(run);
        stat.different_from_chance.(DT_runs{run}).(variable{var}).size_relation.adj_p.size = adjp_size(run);
        end
    end
end
save(stat_path,'stat');
