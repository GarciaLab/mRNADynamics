function plotDorsalResultsLoop(dataType, activityType, varargin)
% valid activites-
% 1. fraction
% 2. timeon
% 3. mRNA
% 4. duration
% 5. max

%%
nc = 12;
paramSearch = [];
R = 1;
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
load([resultsFolder,filesep,DataType,filesep,'dorsalResults.mat'], 'dorsalResults')

nPlots = numel(paramSearch);
cmap = single(summer(nPlots));
cmap2 = single(spring(1));

dorsalResults = dorsalResults{nc-11};
x = dorsalResults.dorsalFluoBins;

y = dorsalResults.fracFluoEmbryo;
ymean = dorsalResults.meanFracFluoEmbryo;
se = dorsalResults.seFracFluoEmbryo;

%%

for plotIndex = 1:nPlots
    
    plotScatter = plotIndex == 1;
    
    param= paramSearch(plotIndex);
    
    plotDorsalActivity(x, y,activityType, nc,...
        dataType, ymean, se, plotScatter,...
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
