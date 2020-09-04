clear
close all

Prefix = '2019-11-26-2xDl_Venus_snaBAC_MCPmCherry_Leica_Zoom45_21uW14uW_01';
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
preProcFolder = liveExperiment.preFolder;
trackFigFolder = [dropboxFolder filesep Prefix '\TrackingFigures\'];
ParticlesFull = getParticlesFull(liveExperiment);

MS2FullParticles = ParticlesFull.FullParticles{1,2};
TFRawParticles = ParticlesFull.RawParticles{1,1};

% For all the nuclei that have final, full MS2 traces, get the xyz position
% info for all raw detections of TF clusters across all time points
nucleiWithMS2Spots = [MS2FullParticles.Nucleus];
clusterNuclei = [TFRawParticles.Nucleus];
numNuclei = numel(nucleiWithMS2Spots);

nucleiWithClusters(numNuclei).nucleusID = [];
nucleiWithClusters(numNuclei).frame = [];
nucleiWithClusters(numNuclei).xPos = [];
nucleiWithClusters(numNuclei).yPos = [];
nucleiWithClusters(numNuclei).zPos = [];
clustersByNucleus(numNuclei).nucleusID = [];
for n = 1:numNuclei
    % Find the clusters in this nucleus
    currNucleus = nucleiWithMS2Spots(n);
    nucleiWithClusters(n).nucleusID = currNucleus;
    clustersByNucleus(n).nucleusID = currNucleus;
    clusterNucleiIndices = find(clusterNuclei == currNucleus);
    
    % Get the frame, xpos, ypos, and zpos of each cluster
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
    
    % Now, sort the positions by frame so that we can access each unique
    % combo of nucleus+frame
    clustersByNucleus(n).frames(nFrames).xyzCoord = [];  %initialize to correct size
    for t = 1:nFrames
        % Figure out whether this nucleus has clusters in this frame
        frameIndices = find(nucleiWithClusters(n).frame == t);
        if ~isempty(frameIndices)
            % Get all the xyz positions for the clusters in this frame
            xyzPos = nan(numel(frameIndices),3);
            for f = 1:numel(frameIndices)
                currIndex = frameIndices(f);
                xyzPos(f,1) = nucleiWithClusters(n).xPos(currIndex);
                xyzPos(f,2) = nucleiWithClusters(n).yPos(currIndex);
                xyzPos(f,3) = nucleiWithClusters(n).zPos(currIndex);
            end
            clustersByNucleus(n).frames(t).xyzCoord = xyzPos;
        end
    end
end