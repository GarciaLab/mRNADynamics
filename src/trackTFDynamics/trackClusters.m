function RawParticles = trackClusters(Prefix,varargin)

% Process user options
% searchRadiusMicrons = 5 by default
[useHistone,searchRadiusMicrons,retrack,displayFigures] = ...
            determinePerformTrackingOptions(varargin);

% Get all the required data for this Prefix
liveExperiment = LiveExperiment(Prefix);

nCh = numel(liveExperiment.spotChannels);   %Only grabbing the spot channels - might cause issues if spot channels aren't the first n channels
pixelSize = liveExperiment.pixelSize_um;    %NL: pixel size is in um
experimentType = liveExperiment.experimentType;
FrameInfo = getFrameInfo(liveExperiment);
Spots = getSpots(liveExperiment);
schnitzCells = getSchnitzcells(liveExperiment);

% [Particles] = track02KalmanTesting(...
%     FrameInfo, Spots, NCh, PixelSize, SearchRadiusMicrons, retrack, displayFigures)
  
tic
disp('Performing intitial particle linking...')
RawParticles = track01ParticleProximity(FrameInfo, Spots, schnitzCells, ...
                nCh, pixelSize, searchRadiusMicrons, useHistone, retrack, ...
                displayFigures);
toc