function plotMS2ClusterDistances(Prefixes)
close all;

% Prefixes = {'2019-11-26-2xDl_Venus_snaBAC_MCPmCherry_Leica_Zoom45_21uW14uW_01'
%             '2020-07-23-2xDl-Ven_hbBAC-MCPmCh_Leica_Zoom45_21uW14uW_10'};
        
% Prefixes = {'2022-03-29-Dl-mNeonGreen-MCP-mCh_snaBAC_settingsTest06_embryo04'
%             '2020-07-23-2xDl-Ven_hbBAC-MCPmCh_Leica_Zoom45_21uW14uW_10'};
countsNorm = {};
edgesNorm = {};
nDetections = [];
for i = 1:numel(Prefixes)
    liveExperiment = LiveExperiment(Prefixes{i});

    nFrames = liveExperiment.nFrames;
    pixelSize_nm = liveExperiment.pixelSize_nm;

    resultsFolder = liveExperiment.resultsFolder;
    clusterResultsFolder = [resultsFolder filesep 'cluster_analysis'];

    load([clusterResultsFolder filesep 'ms2ClusterDistances.mat'],'ms2ClusterDistances');

    % For simplest analysis, just combine data from all nuclei and all frames
    allDistances_um = [];
    numNuclei = numel(ms2ClusterDistances);
    for n = 1:numNuclei
        for f = 1:nFrames
            newDistances = ms2ClusterDistances(n).frames(f).ms2ClusterDist;
            allDistances_um = [allDistances_um, newDistances];
        end
    end

    % The data is in um or px, convert to nm
    % hb data is in px, snaBAC data is in um
%     if i==1
        allDistances_nm = allDistances_um .* 1000;
%     elseif i==2
%         allDistances_px = allDistances_um;
%         allDistances_nm = allDistances_px .* pixelSize_nm;
%     end


    % Find the average number of "interactions" between one MS2 spot and the
    % clusters in that nucleus
    distanceThresh = 180; % in nm

    interactions = allDistances_nm(allDistances_nm <= distanceThresh);
    numInteractions = numel(interactions);
    numInteractionsNormalized = numInteractions/numNuclei/nFrames;  % naive calculation assuming only 1 spot per nucleus

    % Plot distribution of interaction event distances
    histThresh = 5 * distanceThresh;
    histDistances = allDistances_nm(allDistances_nm <= histThresh);

    nBins = ceil(histThresh / pixelSize_nm);
    
    edges = 0:50:900;
    [counts, edges] = histcounts(histDistances,edges);
    countsNorm{i} = counts ./ sum(counts);  %normalized by # of cluster detections to compare across experiments
    edgesNorm{i} = edges;
    nDetections(i) = sum(counts);
end

%% Overlapping histograms
pbocYlw = [234,194,100]/255;
pbocRed = [213,108,85]/255;
h = {};
figOverlay = figure(1);
hold on
h{1} = histogram('BinEdges',edgesNorm{1},'BinCounts',countsNorm{1},'FaceColor', pbocRed);
h{2} = histogram('BinEdges',edgesNorm{2},'BinCounts',countsNorm{2},'FaceColor', pbocYlw);
xline(180)
xlabel('distance of cluster from MS2 spot (nm)')
ylabel('prob')
legend(['snaBAC, n = ' num2str(nDetections(1))], ['hbBAC, n = ' num2str(nDetections(2))]);
hold off

figFolder = 'S:\Meghan\Dropbox\DorsalClustersFigures\cluster_analysis\';
if ~exist(figFolder, 'dir')
    mkdir(figFolder)
end

StandardFigurePBoC(h{1},gca);
exportgraphics(figOverlay,[figFolder filesep 'ms2ClusterDistance_snaVShb_overlay_pboc.png']);
exportgraphics(figOverlay,[figFolder filesep 'ms2ClusterDistance_snaVShb_overlay_pboc.pdf']);


%% Stacked subplots
figSubplotsStack = figure(2);
figSubplotsStack.Position(3:4) = [800,800];
hold on
subplot(2,1,1);
hs{1} = histogram('BinEdges',edgesNorm{1},'BinCounts',countsNorm{1},'FaceColor', pbocRed, 'FaceAlpha', 1);
ylabel('prob')
lgnd1 = legend(['snaBAC, n = ' num2str(nDetections(1))]);
set(lgnd1,'color','none','Box','off','Location','northwest');
StandardFigurePBoC(hs{1},gca);

subplot(2,1,2);
hs{2} = histogram('BinEdges',edgesNorm{2},'BinCounts',countsNorm{2},'FaceColor', pbocYlw, 'FaceAlpha', 1);
xlabel('distance of cluster from MS2 spot (nm)')
ylabel('prob')
lgnd2 = legend(['hbBAC, n = ' num2str(nDetections(2))]);
set(lgnd2,'color','none','Box','off','Location','northwest');
StandardFigurePBoC(hs{2},gca);
hold off

exportgraphics(figSubplotsStack,[figFolder filesep 'ms2ClusterDistance_snaVShb_vertStack_pboc.png']);
exportgraphics(figSubplotsStack,[figFolder filesep 'ms2ClusterDistance_snaVShb_vertStack_pboc.pdf']);

