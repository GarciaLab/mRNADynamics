% script to compare MCP levels across sets
clear
close all

% define Prefix list and plot names
PrefixCell = {...
    '2017-09-14-P2P-MCP-NoNLS-mCherry-doubledosage2',...
    '2017-09-14-P2P-MCP-NoNLS-mCherry-doubledosage3',...
    '2017-12-13-P2P-MCP-NoNLS-mCherry-doubledosage',...
    '2017-09-10-P2P-MCP-NoNLS-mCherry-threeovertwodosage',...
    '2017-09-15-P2P-MCP-NoNLS-mCherry-threeovertwodosage',...
    '2017-09-17-P2P-MCP-NoNLS-mCherry-threeovertwodosage',...
    '2017-12-03-P2P-MCP-NoNLS-mCherry',...
    '2017-12-03-P2P-MCP-NoNLS-mCherry_15',...
    '2019-10-31-F_F_F_dorsal_synthetics_saturation_check',...
    '2019-11-07-F_F_F_MCP-NoNLS_x_HbP2Pv1MS2_01',...
    '2019-11-04-F_F_3_Dlxsna_mcp_check_10umStack'...
    };
% logical vecorrs    
double_vec = [1 1 1 0 0 0 0 0 0 0 0 0];
threetwo_vec = [0 0 0 1 1 1 0 0 0 0 0];
one_vec = [0 0 0 0 0 0 1 1 0 0 0];
fff_vec = [0 0 0 0 0 0 0 0 1 1 0];
ff3_vec = [0 0 0 0 0 0 0 0 0 0 1];
plot_fluo_vec = [1 1 1 1 1 1 1 1 0 1 0];
PlotNameCell = {'2.5x','2.5x','2.5x', '2x','2x','2x', '1x','1x', 'FFF line','FFF line','FF3 line'};
% load data to plot
plot_struct = struct;
plot_times = 1:25; % minutes
ap_range = 25:35;
t_sigma = 1;
n_boots = 100;
for p = 1:numel(PrefixCell)
    Prefix = PrefixCell{p};
    % get paths
    [rawDataPath,ProcPath,DropboxFolder]= readMovieDatabase(Prefix);
    % check to see if MCP datasets exist. If not, generate requisite data
    mcp_path = [DropboxFolder, filesep, Prefix, filesep, 'MCPLevelsData.mat'];
    cp_path = [DropboxFolder, filesep, Prefix, filesep, 'CompiledParticles.mat'];
    if exist(mcp_path)
        % read data
        load(mcp_path, 'MCPLevelsData','MCPLevelsSummary');
    else
        disp(['running MCP analysis for ' Prefix ' ...'])
        MCPLevelAnalysis(Prefix);
        load(mcp_path, 'MCPLevelsData','MCPLevelsSummary');
    end
    plot_struct(p).mcp_data = MCPLevelsSummary;
    plot_struct(p).time = [MCPLevelsData.time];
    %%% incorporate compiledParticles info
    mean_fluo_vec = NaN(size(plot_times));
    mean_fluo_vec_se = NaN(size(plot_times));
    max_fluo = NaN;
    max_fluo_se = NaN;
    fluo_time = NaN;
    if plot_fluo_vec(p) == 1        
        load(cp_path);            
        TimeVec = ElapsedTime(nc14:end) - ElapsedTime(nc14);
        if iscell(AllTracesVector)
            AllTracesVector = AllTracesVector{1};
            CompiledParticles = CompiledParticles{1};
        end
            
        AllTraces14 = AllTracesVector(nc14:end,:);
        act_ft = sum(~isnan(AllTraces14))>5;
        ap_ft = ismember(round(100*[CompiledParticles.MeanAP]),ap_range);
        AllTraces14 = AllTraces14(:,ap_ft&act_ft);
        % bootstrap
        ind_vec = 1:size(AllTraces14,2);       
        if ~isempty(ind_vec)
            mean_fluo_array = NaN(numel(TimeVec),n_boots);
            max_fluo_vec = NaN(1,n_boots);
            for n = 1:n_boots
                s_ids = randsample(ind_vec,numel(ind_vec),true);
                mean_fluo_array(:,n) = nanmean(AllTraces14(:,s_ids),2);
                tr_ft = AllTraces14(:,s_ids);
                max_fluo_vec(n) = prctile(tr_ft(:),95);
            end
            mean_fluo_vec = nanmean(mean_fluo_array,2)';
            mean_fluo_vec_se = nanstd(mean_fluo_array,[],2)';
            max_fluo = nanmean(max_fluo_vec);
            max_fluo_se = nanstd(max_fluo_vec);
            fluo_time = TimeVec;
        else
            warning(['no nc14 trace data for ' Prefix '. Will not plot fluo...'])
            plot_fluo_vec(p) = 0;
        end
    end
    plot_struct(p).mean_fluo_vec = mean_fluo_vec;
    plot_struct(p).mean_fluo_vec_se = mean_fluo_vec_se;
    plot_struct(p).max_fluo = max_fluo;
    plot_struct(p).max_fluo_se = max_fluo_se;
    plot_struct(p).fluo_time = fluo_time;
end
    
%% plot MCP concentration in nculei over time
FigPath = 'E:\Nick\LivemRNA\Dropbox (Personal)\LocalEnrichmentFigures\PipelineOutput\MCP_saturation_checks\';

hm_cm = brewermap(8,'set2');
dbl_color = hm_cm(2,:);
trtwo_color = hm_cm(1,:);
one_color = hm_cm(3,:);
ff3_color = hm_cm(end,:);
plot_color_array = NaN(numel(plot_struct),3);
close all;
med_mcp_fig = figure('Position', [10 10 1024 512]);
hold on
% dobule
lgd_str = {};
lgd_plt = [];
flag = true;
for p = find(double_vec)    
    e2 = errorbar(plot_struct(p).time/60,plot_struct(p).mcp_data.med_mcp_mean,plot_struct(p).mcp_data.med_mcp_se,'Color',dbl_color,'LineWidth',1.5);
    e2.CapSize = 0;
    plot_color_array(p,:) = dbl_color;
    if flag
        lgd_str = [lgd_str{:} {'2x'}];
        lgd_plt = [lgd_plt e2];
        flag = false;
    end
end
% three two
flag = true;
for p = find(threetwo_vec)
    e32 = errorbar(plot_struct(p).time/60,plot_struct(p).mcp_data.med_mcp_mean,plot_struct(p).mcp_data.med_mcp_se,'Color',trtwo_color,'LineWidth',1.5);
    e32.CapSize = 0;
    plot_color_array(p,:) = trtwo_color;
    if flag
        lgd_str = [lgd_str{:} {'3/2x'}];
        lgd_plt = [lgd_plt e32];
        flag = false;
    end    
end
% one
flag = true;
for p = find(one_vec)
    e1 = errorbar(plot_struct(p).time/60,plot_struct(p).mcp_data.med_mcp_mean,plot_struct(p).mcp_data.med_mcp_se,'Color',one_color,'LineWidth',1.5);
    e1.CapSize = 0;
    plot_color_array(p,:) = one_color;
    if flag
        lgd_str = [lgd_str{:} {'1x'}];
        lgd_plt = [lgd_plt e1];
        flag = false;
    end
end

% FFF
flag = true;
for p = find(fff_vec)
    efff = errorbar(plot_struct(p).time/60,plot_struct(p).mcp_data.med_mcp_mean,plot_struct(p).mcp_data.med_mcp_se,'Color','black','LineWidth',1.5);
    efff.CapSize = 0;
    plot_color_array(p,:) = [0 0 0];
    if flag
        lgd_str = [lgd_str{:} {'FFF'}];
        lgd_plt = [lgd_plt efff];
        flag = false;
    end
end

% FF3
flag = true;
for p = find(ff3_vec)
    eff3 = errorbar(plot_struct(p).time/60,plot_struct(p).mcp_data.med_mcp_mean,plot_struct(p).mcp_data.med_mcp_se,'--','Color',ff3_color,'LineWidth',1.5);
    eff3.CapSize = 0;
    plot_color_array(p,:) = ff3_color;
    if flag
        lgd_str = [lgd_str{:} {'FF3'}];
        lgd_plt = [lgd_plt eff3];
        flag = false;
    end
end

xlabel('minutes from start of nc14')
ylabel('median MCP levels (au)')

legend(lgd_plt, lgd_str{:});

box on 
grid on
set(gca,'FontSize',14);

saveas(med_mcp_fig,[FigPath 'mcp_level_comparison.png'])

%% plot cmp level vs max MS2
plot_indices = find(plot_fluo_vec);
plot_time = 5;
max_fluo_fig = figure;
hold on
lgd_str = {};
lgd_plt = [];
plt_ind = [1 4 8 10];
for p = plot_indices
    mf = plot_struct(p).max_fluo;
    mf_se = plot_struct(p).max_fluo_se;
    time = plot_struct(p).time/60;
    [dtmin, ti] = min(abs(time-plot_time));
    if dtmin < 1
%         plt_ind = [plt_ind p];
        mcp = nanmean(plot_struct(p).mcp_data.med_mcp_mean(ti));
        errorbar(mcp,mf,mf_se,'Color','black')
        if ismember(p,plt_ind)
            lgd_plt = [lgd_plt scatter(mcp,mf,'MarkerFaceColor',plot_color_array(p,:),'MarkerEdgeColor','black')];
        else
            scatter(mcp,mf,'MarkerFaceColor',plot_color_array(p,:),'MarkerEdgeColor','black');
        end
    end
end
grid on
box on
ylim([200 1000])
xlabel('nuclear MCP concentration (au)')
ylabel('hbP2Pv1 95th prectentile (au)')
set(gca,'Fontsize',16);
legend(lgd_plt,PlotNameCell{plt_ind},'Location','southeast')
saveas(max_fluo_fig,[FigPath 'mcp_level_vs_hbP2P.png'])