function plotDorsalResultsLoop(dataType,...
    activityType, paramRange, modelType, varargin)
% valid activites-
% 1. fraction
% 2. timeon
% 3. mRNA
% 4. duration
% 5. max
%
% paramSearch %list of parameter values to plot. keep small (<10)
%
%modelType- string containing a kind of model
% 1.hill
% 2.simpleWithPol
% 3. mwcNoPol
%
arguments
    dataType char
    activityType char
    paramRange double
    modelType char
end
arguments(Repeating)
    varargin
end


%%
nc = 12;
R = 1; %the rate/amplitude parameter of the models in fitDorsalActivity
xRange = [0 4000];
legendVisible = 'off';
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

nPlots = numel(paramRange);
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
    
    paramValue = paramRange(plotIndex);
    
    plotDorsalActivity(dorsalFluoBins, dorsalActivity, activityType, nc,...
        dataType, dorsalActivity_mean, dorsalActivity_SE, shouldPlotScatter,...
        'modelType', modelType,...
        'fix1', R, 'fix4', 0, 'fix5', paramValue);
    
    if plotIndex == 1
        ax1 = plotInLoop(plotIndex, cmap2,...
            'xRange', xRange, 'legendVisible', legendVisible);
    else
        plotInLoop(plotIndex, cmap, 'xRange',...
            xRange, 'ax1', ax1, 'legendVisible', legendVisible);
    end
    
end


end
