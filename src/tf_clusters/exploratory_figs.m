function exploratory_figs(prefix)
close all;

%% Get prefix info
liveExperiment = LiveExperiment(prefix);

% folders
dropboxFolder = liveExperiment.userResultsFolder;
outputFolder = [dropboxFolder, filesep, prefix];
clusterResultsFolder = [outputFolder, filesep, 'cluster_analysis'];

% channel info
inputChannels = liveExperiment.inputChannels;
spotChannels = liveExperiment.spotChannels;
nSpotCh = numel(spotChannels);
clusterChannel = intersect(inputChannels, spotChannels);
ms2Channel = setxor(inputChannels, spotChannels);

%% Load in cluster data
[Ellipses, schnitzcells] = loadNucleiSegmentationAndLineages(outputFolder, prefix);

% [Spots, ~] = loadSpotsAndCreateSpotFilter(dropboxFolder, prefix, nSpotCh);
% SpotsMS2 = Spots{1,ms2Channel};
% SpotsClusters = Spots{1,clusterChannel};

Clusters = loadClusters(outputFolder);

[Particles, ~] = getParticles(liveExperiment);
Particles = Particles{1,ms2Channel};


%% Visualize initial segmentation results
nucleusID = 5;
frame = 60;
ch = 2;

movieFrame = getMovieFrame(liveExperiment, frame, ch);
frameMax = max(movieFrame,3);

clustersX = [Clusters(nucleusID).ClustersByFrame(frame).ClusterFits.xPos];
clustersY = [Clusters(nucleusID).ClustersByFrame(frame).ClusterFits.yPos];

particleNucIndex = find([Particles.Nucleus]==nucleusID);
particleFrameIndex = find(Particles(particleNucIndex).Frame==frame);
ms2X = Particles(particleNucIndex).xPos(particleFrameIndex);
ms2Y = Particles(particleNucIndex).yPos(particleFrameIndex);

figure
imshow(frameMax, []);
hold on
scatter(clustersX, clustersY, 50, 'green','filled');
scatter(ms2X, ms2Y, 50, 'red','filled');
hold off




%% Final data from ms2ClusterDistances
% load([clusterResultsFolder, filesep, 'ms2ClusterDistances.mat'],...
% 'ms2ClusterDistances');
% % ms2ClusterDistances = ms2ClusterDistances_sna;
% 
% clustersXPos = ms2ClusterDistances(nucleusID).frames(frame).xyzCoordClusters(:,1);
% clustersYPos = ms2ClusterDistances(nucleusID).frames(frame).xyzCoordClusters(:,2);
% clustersZPos = ms2ClusterDistances(nucleusID).frames(frame).xyzCoordClusters(:,3);
% 
% ms2XPos = ms2ClusterDistances(nucleusID).frames(frame).xyzCoordMS2(1);
% ms2YPos = ms2ClusterDistances(nucleusID).frames(frame).xyzCoordMS2(2);
% ms2ZPos = ms2ClusterDistances(nucleusID).frames(frame).xyzCoordMS2(3);
% 
% figure(2)
% scatter3(clustersXPos,clustersYPos,clustersZPos, 50, 'black')
% hold on
% scatter3(ms2XPos,ms2YPos,ms2ZPos, 50, 'red')
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% hold off
