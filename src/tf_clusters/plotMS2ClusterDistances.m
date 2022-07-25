function plotMS2ClusterDistances(dataTypes) %,clusterTF)
%
% DESCRIPTION
%
%
% INPUT ARGUMENTS
% dataTypes: This is a cell array of char variable(s) which are the exact 
%            name(s) of the tab(s) in DataStatus.xlsx that you wish to 
%            analyze. Single dataType can also be passed as a string (char
%            array).
%            NB: This function expects dataStatus tabs to be named in a
%                specific way: the transcription factor, then the 
%                MS2-tagged gene, separated by an underscore. 
%            e.g. dataTypes = {'Dl-mNG_snaBAC', 'Dl-mNG_hbP2P', 'DleGFP_twiPEe'}
%
% clusterTF: A string specifying which TF has the clusters, used to determine
%            which genes are target vs. non-target
%
%
% OPTIONS
%
%
% OUTPUT
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 2020/08/12
% Last Updated: 2022/05/24
%

% Close any open figures so we aren't drawing on top of other data
close all;

%% Parse input arguments
% If the user passed a single dataType as a string, convert to cell array 
% for downstream compatibility
if ischar(dataTypes)
    dataTypeList = {dataTypes};
else
    dataTypeList = dataTypes;
end

% Get the prefixes for each dataType, and keep them organized by dataType
% compiledDistances
distancesByPrefix = struct('DataType',[],'Prefix',[],'embryoID',[], ...
                           'all_distances_nm',[]);
prefix_counter = 1;
for i = 1:length(dataTypeList)
    currDataType = dataTypeList{i};

    currPrefixes = getProjectPrefixes(currDataType, 'customApproved',...
                                        'ReadyForMS2ClusterDistancePlots');
    for j = 1:numel(currPrefixes)
        thisPrefix = currPrefixes{j};
        
        % Add the datatype and prefix names to main structure
        distancesByPrefix(prefix_counter).DataType = currDataType;
        distancesByPrefix(prefix_counter).Prefix = thisPrefix;
        
        % We want to combine the datatype + embryoID to make identifiers 
        % for each experiment, so we need to extract the embryoID
        strStartIndex = strfind(thisPrefix, 'embryo');
        if numel(strStartIndex) == 1
            embryoID = thisPrefix(strStartIndex:strStartIndex+7);
            distancesByPrefix(prefix_counter).embryoID = embryoID;
        else
            error('The string ''embryo'' appears more than once. Can''t handle this case.');
        end
        
        prefix_counter = prefix_counter+1;
    end
end
numPrefixes = numel(distancesByPrefix);

%% Set parameters for generic histogram plots
% Our detection limit is roughly 200nm (anything below that is in the
% same space, as far as we can tell), so let's use that to set some
% reasonable parameters on the data we plot
interactionThresh = 200; % in nm

%     interactions = distances_nm(distances_nm <= distanceThresh);
%     numInteractions = numel(interactions);
%     numInteractionsNormalized = numInteractions/numNuclei/nFrames;  % naive calculation assuming only 1 spot per nucleus

% Standard radial enrichment plots in the literature usually only go out to
% 1um, so we'll start with that.
maxDistanceForHist = 1000; % in nm

%     nBins = ceil(maxDistanceForHist / pixelSize_nm);

% I should do this in a more principled way, using xy- & z-step sizes or
% interpolating
% For now, use common step size in other enrichment papers, 0.05um
binStepSize = 50;
% inputEdges = 0:binStepSize:maxDistanceForHist;
% inputEdges = 0:binStepSize:3000;
inputEdges = 0:binStepSize:maxDistanceForHist;

% Need to correct for differing volumes in each bin (each bin is actually a
% a spherical shell because we're doing this in 3D
binVolumes3D = nan(1, numel(inputEdges)-1);
V = @(r) (4/3)*pi*(r^3);

for i = 1:numel(inputEdges)-1
    binVolumes3D(i) = V(inputEdges(i+1)) - V(inputEdges(i));
end

binVolumes3D = binVolumes3D ./ V(inputEdges(end)); %just rescaling so we're not working with tiny, tiny counts

 
%% Compile stats for each prefix so we can plot them seperately
for p = 1:numPrefixes
    % Get relevant info for each prefix
    liveExperiment = LiveExperiment(distancesByPrefix(p).Prefix);
    nFrames = liveExperiment.nFrames;

    clusterResultsFolder = [liveExperiment.resultsFolder, filesep, ...
                            'cluster_analysis'];
    load([clusterResultsFolder, filesep, 'ms2ClusterDistances.mat'],...
         'ms2ClusterDistances');

    % For simplest analysis, just combine data from all nuclei and all frames
    distances_um = [];
    numNuclei = numel(ms2ClusterDistances);
    for n = 1:numNuclei
        for f = 1:nFrames
            newDistances = ms2ClusterDistances(n).frames(f).ms2ClusterDist;
            distances_um = [distances_um, newDistances];
        end
    end
    
    % The data is in um, convert to nm
    distances_nm = distances_um .* 1000;

%     distances_nm = distances_um * pixelSize_nm; %these are actually in pixels
    
    % Bin data to make histogram plot
    histDistances = distances_nm(distances_nm <= maxDistanceForHist);
%     histDistances = distances_nm;
    [histCountsRaw, histEdges] = histcounts(histDistances,inputEdges);
    
    histCountsPerVoxel = histCountsRaw ./ binVolumes3D;
    
    histCountsPerVoxelNorm = histCountsPerVoxel ./ sum(histCountsRaw);  %normalized by # of cluster detections to compare across experiment

    % add all data to main structure
    distancesByPrefix(p).all_distances_nm = distances_nm;
    distancesByPrefix(p).histEdges = histEdges;
    distancesByPrefix(p).histCountsRaw = histCountsRaw;
    distancesByPrefix(p).histCountsPerVoxel = histCountsPerVoxel;
    distancesByPrefix(p).histCountsPerVoxelNorm = histCountsPerVoxelNorm;
end

%% Combine all prefixes of the same DataType to plot together
% TODO
distancesByDataType = struct('DataType',dataTypes,'distances_all_nm',[]);
% 
for p = 1:numel(distancesByPrefix)
    
    switch distancesByPrefix(p).DataType
        case dataTypes{1}
            distancesByDataType(1).distances_all_nm = [distancesByDataType(1).distances_all_nm, distancesByPrefix(p).all_distances_nm];
        case dataTypes{2}
            distancesByDataType(2).distances_all_nm = [distancesByDataType(2).distances_all_nm, distancesByPrefix(p).all_distances_nm];
    end
    
end

for d = 1:numel(distancesByDataType)
    histDistances = distancesByDataType(d).distances_all_nm(distancesByDataType(d).distances_all_nm <= maxDistanceForHist);
    [histCountsRaw, histEdges] = histcounts(histDistances,inputEdges);

    histCountsPerVoxel = histCountsRaw ./ binVolumes3D;

    histCountsPerVoxelNorm = histCountsPerVoxel ./ sum(histCountsRaw);  %normalized by # of cluster detections to compare across experiment

    % add all data to main structure
    distancesByDataType(d).histEdges = histEdges;
    distancesByDataType(d).histCountsRaw = histCountsRaw;
    distancesByDataType(d).histCountsPerVoxel = histCountsPerVoxel;
    distancesByDataType(d).histCountsPerVoxelNorm = histCountsPerVoxelNorm;
end

%% Generic plotting elements
todaysDate = datestr(now, 'yyyy-mm-dd');
figFolder = ['S:\Meghan\Dropbox\DorsalClustersFigures\cluster_analysis\',...
             todaysDate];
if ~exist(figFolder, 'dir')
    mkdir(figFolder)
end

pbocYlw = [234,194,100]/255;
pbocRed = [213,108,85]/255;

targetGeneColor = pbocRed;
nontargetGeneColor = pbocYlw;

%% Individual histograms 
for p = 1:numPrefixes
    figIndivHist = figure(p);
    figIndivHist.Position = [500,500,1000,330];
    
    if contains(distancesByPrefix(p).DataType, 'sna', 'IgnoreCase', true)
        histColor = targetGeneColor;
    elseif contains(distancesByPrefix(p).DataType, 'hb', 'IgnoreCase', true)
        histColor = nontargetGeneColor;
    else
        error('Not a supported gene.');
    end
    expName = [distancesByPrefix(p).DataType,'_', distancesByPrefix(p).embryoID];
    nDetections = sum([distancesByPrefix(p).histCountsRaw]);
    
    hold on
    hst = histogram('BinEdges',distancesByPrefix(p).histEdges,...
                      'BinCounts',distancesByPrefix(p).histCountsPerVoxelNorm,...
                      'FaceColor', histColor, 'FaceAlpha', 1);
    ylim([0 15])
    ylabel({'number of cluster detections';'per voxel, normalized'})
    xlabel('distance from MS2 locus (nm)')
    lgnd = legend([expName, ', n = ' num2str(nDetections)],'Interpreter','none');
    set(lgnd,'color','none','Box','off','Location','northwest');
    StandardFigurePBoC(hst,gca);
    hold off;
    
    % Save figure

    filename = ['ms2ClusterDistance_', expName]; 
    exportgraphics(figIndivHist,[figFolder, filesep, filename, '.png']);
    exportgraphics(figIndivHist,[figFolder, filesep, filename, '.pdf']);
end

%% Combine all embryos by DataType histograms
for d = 1:numel(distancesByDataType)
    figIndivHist = figure(p+d);
    figIndivHist.Position = [500,500,1000,330];
    
    if contains(distancesByDataType(d).DataType, 'sna', 'IgnoreCase', true)
        histColor = targetGeneColor;
    elseif contains(distancesByDataType(d).DataType, 'hb', 'IgnoreCase', true)
        histColor = nontargetGeneColor;
    else
        error('Not a supported gene.');
    end
    expName = [distancesByDataType(d).DataType,'_allembryos'];
    nDetections = sum([distancesByDataType(d).histCountsRaw]);
    
    hold on
    hst = histogram('BinEdges',distancesByDataType(d).histEdges,...
                      'BinCounts',distancesByDataType(d).histCountsPerVoxelNorm,...
                      'FaceColor', histColor, 'FaceAlpha', 1);
    ylim([0 15])
    ylabel({'number of cluster detections';'per voxel, normalized'})
    xlabel('distance from MS2 locus (nm)')
    lgnd = legend([expName, ', n = ' num2str(nDetections)],'Interpreter','none');
    set(lgnd,'color','none','Box','off','Location','northwest');
    StandardFigurePBoC(hst,gca);
    hold off;
    
    % Save figure

    filename = ['ms2ClusterDistance_', expName]; 
    exportgraphics(figIndivHist,[figFolder, filesep, filename, '.png']);
    exportgraphics(figIndivHist,[figFolder, filesep, filename, '.pdf']);
end

%% Stacked subplots
% figSubplotsStack = figure(2);
% figSubplotsStack.Position(3:4) = [800,800];
% hold on
% subplot(2,1,1);
% hs{1} = histogram('BinEdges',edgesNorm{1},'BinCounts',countsNorm{1},'FaceColor', pbocRed, 'FaceAlpha', 1);
% ylabel('prob')
% lgnd1 = legend(['snaBAC, n = ' num2str(nDetections(1))]);
% set(lgnd1,'color','none','Box','off','Location','northwest');
% StandardFigurePBoC(hs{1},gca);
% 
% subplot(2,1,2);
% hs{2} = histogram('BinEdges',edgesNorm{2},'BinCounts',countsNorm{2},'FaceColor', pbocYlw, 'FaceAlpha', 1);
% xlabel('distance of cluster from MS2 spot (nm)')
% ylabel('prob')
% lgnd2 = legend(['hbBAC, n = ' num2str(nDetections(2))]);
% set(lgnd2,'color','none','Box','off','Location','northwest');
% StandardFigurePBoC(hs{2},gca);
% hold off
% 
% exportgraphics(figSubplotsStack,[figFolder filesep 'ms2ClusterDistance_snaVShb_vertStack_pboc.png']);
% exportgraphics(figSubplotsStack,[figFolder filesep 'ms2ClusterDistance_snaVShb_vertStack_pboc.pdf']);

