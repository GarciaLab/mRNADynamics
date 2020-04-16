function plotDorsalResultsLoop(dataType,...
    activityType, paramSearch)
% valid activites-
% 1. fraction
% 2. timeon
% 3. mRNA
% 4. duration
% 5. max
%
% paramSearch %list of parameter values to plot. keep small (<10)

arguments
    
    dataType char
    activityType char
    paramSearch double
end


%%
nc = 12;
R = 1; %the rate/amplitude parameter of the models in fitDorsalActivity
xRange = [0 4000];
legendVisible = 'off';
modelType = 'hill';
%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

%%
[~, resultsFolder, ~] = getDorsalPrefixes(dataType);
load([resultsFolder,filesep,dataType,filesep,'dorsalResults.mat'], 'dorsalResults')

nPlots = numel(paramSearch);
cmap = single(summer(nPlots));
cmap2 = single(spring(1));

dorsalResults = dorsalResults{nc-11};
dorsalFluoBins = dorsalResults.dorsalFluoBins;

dorsalActivity = dorsalResults.fracFluoEmbryo;
dorsalActivity_mean = dorsalResults.meanFracFluoEmbryo;
dorsalActivity_SE = dorsalResults.seFracFluoEmbryo;

%%

for plotIndex = 1:nPlots
    
    %the data gets a scatter plot and
    %fits get a line plot
    shouldPlotScatter= plotIndex == 1;
    
    param = paramSearch(plotIndex);
    
    plotDorsalActivity(dorsalFluoBins, dorsalActivity, activityType, nc,...
        dataType, dorsalActivity_mean, dorsalActivity_SE, shouldPlotScatter,...
        'modelType', modelType,...
        'fix1', R, 'fix4', 0, 'fix5', param);
    
    if plotIndex == 1
        ax1 = plotInLoop(plotIndex, cmap2,...
            'xRange', xRange, 'legendVisible', legendVisible);
    else
        plotInLoop(plotIndex, cmap, 'xRange',...
            xRange, 'ax1', ax1, 'legendVisible', legendVisible);
    end
    
end


end
