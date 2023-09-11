%% Calculate relation outcome measures and size of the map 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% general info

stat_path = 'stat_NumerosityAll_TimingAll_minVEtime=0.2_minVEnum=0.3_whichCombi=4_topo_measure=mean.mat';

minAmountComparison = 1;
load(stat_path)
DT_run_1 = 'NumerosityAll';
DT_run_2 = 'TimingAll';
DT_runs = {DT_run_1, DT_run_2};

if contains(DT_runs{1},"Numerosity")
    par_1 ={'x0'};
else
    par_1={'x0','y0'};
end

if contains(DT_runs{2},"Numerosity")
    par_2 ={'x0'};
else
    par_2={'x0','y0'};
end

%% analysis per participant
ps_size = [];  ps_hemi = []; ps_subj = [];

for run = 1:length(DT_runs)
    for p1 = 1:length(par_1)
        for p2 = 1:length(par_2)

            not_nan = ~isnan(stat.data.(DT_runs{run}).topo_std.(par_1{p1}).(par_2{p2}));
            size = nan(length(stat.data.(DT_runs{run}).topo_std.(par_1{p1}).(par_2{p2})),1);

            for datapoint = 1:length(stat.data.(DT_runs{run}).topo_std.(par_1{p1}).(par_2{p2}))
                eval(['size(datapoint) = length([all_combis.shared_coord_ids{all_combis.subject == stat.data.(DT_runs{run}).subj(datapoint) & all_combis.hemi == stat.data.(DT_runs{run}).hemi(datapoint) & all_combis.Map', num2str(run) ,'== stat.data.(DT_runs{run}).map(datapoint) & ~contains(all_combis.Map', num2str(3-run) ',"All")}]);'])
            end

            maps = unique(stat.data.(DT_runs{run}).map);
            for map = 1:length(maps)
                if sum(~isnan(stat.data.(DT_runs{run}).topo_std.x0.x0(stat.data.(DT_runs{run}).map==maps{map}))) <= minAmountComparison
                    stat.data.(DT_runs{run}).topo_std.(par_1{p1}).(par_2{p2})(stat.data.(DT_runs{run}).map==maps{map}) = NaN;
                end
            end


            covariates = '{size';
            covariate_names = {'size'};

            covariates = [covariates,', stat.data.(DT_runs{run}).hemi'];
            covariate_names{end+1} = 'hemi';


            covariates = [covariates,', stat.data.(DT_runs{run}).subj'];
            covariate_names{end+1} = 'subj';
            covariates = [covariates, '}'];

            [~, tbl, stats, ~] = anovan(stat.data.(DT_runs{run}).topo_std.(par_1{p1}).(par_2{p2}), eval(covariates),'Continuous',1,'varnames',covariate_names);
            stat.different_from_chance.(DT_runs{run}).topo_per_subj.(par_1{p1}).(par_2{p2}).size_relation.tbl = tbl;
            stat.different_from_chance.(DT_runs{run}).topo_per_subj.(par_1{p1}).(par_2{p2}).size_relation.stats = stats;

            ps_size = [ps_size stat.different_from_chance.(DT_runs{run}).topo_per_subj.(par_1{p1}).(par_2{p2}).size_relation.tbl{2,7}];
            ps_hemi = [ps_hemi stat.different_from_chance.(DT_runs{run}).topo_per_subj.(par_1{p1}).(par_2{p2}).size_relation.tbl{3,7}];
            ps_subj = [ps_subj stat.different_from_chance.(DT_runs{run}).topo_per_subj.(par_1{p1}).(par_2{p2}).size_relation.tbl{4,7}];
        end

    end

end

[~,~,~,adj_p_size] = fdr_bh(ps_size);
[~,~,~,adj_p_hemi] = fdr_bh(ps_hemi);
[~,~,~,adj_p_subj] = fdr_bh(ps_subj);

n = 0;
for run = 1:length(DT_runs)
    for p1 = 1:length(par_1)
        for p2 = 1:length(par_2)
            n = n+1;
            stat.different_from_chance.(DT_runs{run}).topo_per_subj.(par_1{p1}).(par_2{p2}).size_relation.adj_p.size = adj_p_size(n);
            stat.different_from_chance.(DT_runs{run}).topo_per_subj.(par_1{p1}).(par_2{p2}).size_relation.adj_p.hemi = adj_p_hemi(n);
            stat.different_from_chance.(DT_runs{run}).topo_per_subj.(par_1{p1}).(par_2{p2}).size_relation.adj_p.subj =adj_p_subj(n);
        end
    end
end


%% analysis BF over the entire group
BF.NumerosityAll.x0.x0 = [0.3358220,3.7892062,0.8212109,0.3130050,0.4090667,4.2467176,1.6032793, 1.6517540];
BF.NumerosityAll.x0.y0 = [0.4700467,0.6256568,0.7962206,0.2606131,0.4692343,282.1445297,0.9649128,2.2703932];
BF.TimingAll.x0.x0 = [0.3358220,3.7892062,0.6948372,0.2707176,0.2695510,1.4263719,1.6032793,0.3305964,0.6731898];
BF.TimingAll.x0.y0 = [0.4700467,0.6256568,0.6788992,0.2612219,0.2614731,44.2828100,0.9649128,0.3422824,1.8632803];

BF_rois.NumerosityAll = {'NFI','NFS','NLO','NPCI','NPCM','NPCS','NPO','NTO'};
BF_rois.TimingAll = {'TFI','TFS','TLO','TPCI','TPCM','TPCS','TPO','TTOA','TTOP'};

ps=[];
for run = 1:length(DT_runs)
    for p1 = 1:length(par_1)
        for p2 = 1:length(par_2)

            maps = BF_rois.(DT_runs{run});
            for map = 1:length(maps)
                if sum(~isnan(stat.data.(DT_runs{run}).topo_std.x0.x0(stat.data.(DT_runs{run}).map==maps{map}))) <= minAmountComparison
                    size(stat.data.(DT_runs{run}).map==maps{map}) = NaN;
                end
                mean_size(map) = mean(size(stat.data.(DT_runs{run}).map==maps{map}),'all','omitnan');
            end

            [r, p] = corr(mean_size',BF.(DT_runs{run}).(par_1{p1}).(par_2{p2})','type','Spearman');

            stat.different_from_chance.(DT_runs{run}).topo_group_BF.(par_1{p1}).(par_2{p2}).size_relation.r = r;
            stat.different_from_chance.(DT_runs{run}).topo_group_BF.(par_1{p1}).(par_2{p2}).size_relation.p = p;

            ps = [ps p];
        end
    end

end

[~,~,~,adj_ps] = fdr_bh(ps);
nn = 0;
for run = 1:length(DT_runs)
    for p1 = 1:length(par_1)
        for p2 = 1:length(par_2)
            nn = nn +1;
            stat.different_from_chance.(DT_runs{run}).topo_group_BF.(par_1{p1}).(par_2{p2}).size_relation.adj_p = adj_ps(nn);
        end
    end
end

save(stat_path,'stat');
