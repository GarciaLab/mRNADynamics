% script to compare MCP levels across sets
clear
close all

% define Prefix list and plot names
PrefixCell = {'2017-09-14-P2P-MCP-NoNLS-mCherry-doubledosage2','2017-09-14-P2P-MCP-NoNLS-mCherry-doubledosage3','2017-12-13-P2P-MCP-NoNLS-mCherry-doubledosage','2017-09-10-P2P-MCP-NoNLS-mCherry-threeovertwodosage','2017-09-15-P2P-MCP-NoNLS-mCherry-threeovertwodosage','2017-09-17-P2P-MCP-NoNLS-mCherry-threeovertwodosage',...
    '2017-12-03-P2P-MCP-NoNLS-mCherry','2017-12-03-P2P-MCP-NoNLS-mCherry_15','2019-10-31-F_F_F_dorsal_synthetics_saturation_check'};
    
PlotNameCell = {'2.5x (1)','2.5x (2)','2.5x (3)', '2x (1)','2x (2)','2x (3)', '1x (1)','1x (2)', 'FFF line'};
% load data to plot
plot_struct = struct;
for p = 1:numel(PrefixCell)
    Prefix = PrefixCell{p};
    % get paths
    [rawDataPath,ProcPath,DropboxFolder]= readMovieDatabase(Prefix);
    % check to see if MCP datasets exist. If not, generate requisite data
    f_path = [DropboxFolder, filesep, Prefix, filesep, 'MCPLevelsData.mat'];
    if exist(f_path)
        % read data
        load(f_path, 'MCPLevelsData','MCPLevelsSummary');
    else
        disp(['running MCP analysis for ' Prefix ' ...'])
        MCPLevelAnalysis(Prefix);
        load(f_path, 'MCPLevelsData','MCPLevelsSummary');
    end
    plot_struct(p).mcp_data = MCPLevelsSummary;
    plot_struct(p).time = [MCPLevelsData.time];
end
    
%% make plots

med_mcp_fig = figure('Position', [10 10 1024 512]);
hold on
for p = 1:numel(PrefixCell)
    e = errorbar(plot_struct(p).time/60,plot_struct(p).mcp_data.med_mcp_mean,plot_struct(p).mcp_data.med_mcp_se,'LineWidth',1.5);
    e.CapSize = 0;
end
xlabel('minutes from start of nc14')
ylabel('median MCP levels (au)')
if ~isempty(PlotNameCell)
    legend(PlotNameCell{:});
else
    legend(PrefixCell{:});
end
box on 
grid on

set(gca,'FontSize',14);


min_mcp_fig = figure;
hold on
for p = 1:numel(PrefixCell)
    e = errorbar(plot_struct(p).time/60,plot_struct(p).mcp_data.min_mcp_mean,plot_struct(p).mcp_data.min_mcp_se,'LineWidth',1.5);
    e.CapSize = 0;
end
xlabel('minutes from start of nc14')
ylabel('minimum MCP levels (au)')
if ~isempty(PlotNameCell)
    legend(PlotNameCell{:});
else
    legend(PrefixCell{:});
end
box on 
grid on

set(gca,'FontSize',14);