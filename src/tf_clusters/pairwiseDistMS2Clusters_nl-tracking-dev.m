function ms2ClusterDistances = pairwiseDistMS2Clusters_nl-tracking-dev(Prefix)
%
% DESCRIPTION
% This function takes segmented MS2 spots and TF clusters (using the
% Spots.m structure) and calculates the pairwise distances between the 
% MS2 spot in a nucleus and all clusters detected in the same nucleus.
%
% Note: This function currently only works for 1spot data
%
%
% INPUT ARGUMENTS
% Prefix: prefix for this experiment
% 
% OPTIONS
% N/A
%
% OUTPUT
% ms2ClusterDistances: A structure containing the pairwise distances  
%                      between an MS2 spot and all TF clusters detected,
%                      for each nucleus and time frame.
%
% ms2ClusterDistances.mat: MATLAB file containing the ms2ClusterDistances
%                          structure, saved to a subdirectory, entitled
%                          'clusterAnalysis', within your Results folder
%                          (aka, your 'DropboxFolder')
%
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 08/12/2020
% Last Updated: 11/08/2021
%

liveExperiment = LiveExperiment(Prefix);

nCh = numel(liveExperiment.spotChannels);
channels = liveExperiment.Channels;
channelNames = cell(1,nCh);

xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
zDim = liveExperiment.zDim;
nFrames = liveExperiment.nFrames;

Spots = getSpots(liveExperiment);
if ~iscell(Spots)% NL: added for backwards compatibility
  Spots = {Spots};
end

dropboxFolder = liveExperiment.userResultsFolder;
resultsFolder = liveExperiment.resultsFolder;
writeFolder = [resultsFolder filesep 'clusterAnalysis'];
if ~exist(writeFolder,'dir')
    mkdir(writeFolder)
end
preProcFolder = liveExperiment.preFolder;
trackFigFolder = [dropboxFolder filesep Prefix '\TrackingFigures\'];
ParticlesFull = getParticlesFull(liveExperiment);

MS2FullParticles = ParticlesFull.FullParticles{1,2};
TFRawParticles = ParticlesFull.RawParticles{1,1};

%% For all the nuclei that have final, full MS2 traces, get the xyz position
% info for all raw detections of TF clusters across all time points
nucleiWithMS2Spots = [MS2FullParticles.Nucleus];
clusterNuclei = [TFRawParticles.Nucleus];
numNuclei = numel(nucleiWithMS2Spots);

% pre-allocate space in the structures
nucleiWithClusters(numNuclei).nucleusID = [];
nucleiWithClusters(numNuclei).frame = [];
nucleiWithClusters(numNuclei).xPos = [];
nucleiWithClusters(numNuclei).yPos = [];
nucleiWithClusters(numNuclei).zPos = [];
clustersByNucleus(numNuclei).nucleusID = [];
ms2ByNucleus(numNuclei).nucleusID = [];

% Loop over all nuclei with MS2 spots
for n = 1:numNuclei
    currNucleus = nucleiWithMS2Spots(n);
    nucleiWithClusters(n).nucleusID = currNucleus;
    clustersByNucleus(n).nucleusID = currNucleus;
    ms2ByNucleus(n).nucleusID = currNucleus;
    
    % Find the clusters in this nucleus
    clusterNucleiIndices = find(clusterNuclei == currNucleus);
    
    % Get the frame, xpos, ypos, and zpos of each cluster in this nucleus
    for i = 1:numel(clusterNucleiIndices)
        currIndex = clusterNucleiIndices(i);
        
        currFrames = TFRawParticles(currIndex).Frame;
        currXPos = TFRawParticles(currIndex).xPos;
        currYPos = TFRawParticles(currIndex).yPos;
        currZPos = TFRawParticles(currIndex).zPos;
        
        nucleiWithClusters(n).frame = [nucleiWithClusters(n).frame, currFrames];
        nucleiWithClusters(n).xPos = [nucleiWithClusters(n).xPos, currXPos];
        nucleiWithClusters(n).yPos = [nucleiWithClusters(n).yPos, currYPos];
        nucleiWithClusters(n).zPos = [nucleiWithClusters(n).zPos, currZPos];
    end
    
    % Now, sort the xyz positions by frame so that we can compare to the
    % corresponding MS2 spot
    clustersByNucleus(n).frames(nFrames).xyzCoord = [];  %initialize to correct size
    for t = 1:nFrames
        % Figure out whether this nucleus has clusters in this frame
        frameIndices = find(nucleiWithClusters(n).frame == t);
        if ~isempty(frameIndices)
            % Get all the xyz positions for the clusters in this frame
            xyzPos = nan(numel(frameIndices),3);
            for f = 1:numel(frameIndices)
                currIndex = frameIndices(f);
                % putting the xyz coords in a structure that will be easy 
                % to pass directly to pdist2 later on in this script
                xyzPos(f,1) = nucleiWithClusters(n).xPos(currIndex);
                xyzPos(f,2) = nucleiWithClusters(n).yPos(currIndex);
                xyzPos(f,3) = nucleiWithClusters(n).zPos(currIndex);
            end
            
            % Remove duplicate spot detections to minimize redundancy
            % (not sure why these aren't getting filtered out in
            % segmentSpots ... MT: TODO)
            uniqueXYZPos = unique(xyzPos,'stable','rows');
            
            % Store the coordinates by nucleus, frame
            clustersByNucleus(n).frames(t).xyzCoord = uniqueXYZPos;
        end
    end
    
    % Get the MS2 spot coordinates into a similarly useable structure
    ms2ByNucleus(n).frames(nFrames).xyzCoord = []; %initialize to correct size
    currMS2Frames = MS2FullParticles(n).Frame;
    for j = 1:numel(currMS2Frames)
        currIndex = currMS2Frames(j);
        
        % putting the xyz coords in a structure that will be easy 
        % to pass directly to pdist2 later on in this script
        ms2Pos(1,1) = MS2FullParticles(n).xPos(j);
        ms2Pos(1,2) = MS2FullParticles(n).yPos(j);
        ms2Pos(1,3) = MS2FullParticles(n).zPos(j);
        ms2ByNucleus(n).frames(currIndex).xyzCoord = ms2Pos;
    end
end


%% Calculate distances between the MS2 spot and all clusters in each frame 
% and each nucleus
ms2ClusterDistances(numNuclei).nucleusID = [];  %pre-allocated for speed
for n = 1:numNuclei
    currNucleus = nucleiWithMS2Spots(n);
    ms2ClusterDistances(n).nucleusID = currNucleus;
    
    for f = 1:nFrames
        ms2ClusterDistances(n).frames(f).frameID = f;

        ms2Coord = ms2ByNucleus(n).frames(f).xyzCoord;
        clusterCoords = clustersByNucleus(n).frames(f).xyzCoord;
        
        % Pairwise distances
        if ~isempty(ms2Coord) && ~isempty(clusterCoords)
            ms2ClusterDist = pdist2(ms2Coord, clusterCoords);
            ms2ClusterDistances(n).frames(f).ms2ClusterDist = ms2ClusterDist;
        else
            ms2ClusterDistances(n).frames(f).ms2ClusterDist = [];   % make sure it exists even if there's nothing there
        end
    end
end

save([writeFolder filesep 'ms2ClusterDistances.mat'],'ms2ClusterDistances');