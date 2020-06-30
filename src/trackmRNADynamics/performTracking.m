function Particles = performTracking(Prefix,varargin)

% Process user options
% searchRadiusMicrons = 5; by default
[useHistone,searchRadiusMicrons,retrack,displayFigures] = ...
            determinePerformTrackingOptions(varargin);

% Get all the required data for this Prefix
liveExperiment = LiveExperiment(Prefix);

nCh = numel(liveExperiment.spotChannels);   %Only grabbing the spot channels - might cause issues if spot channels aren't the first n channels
pixelSize = liveExperiment.pixelSize_um;    %NL says pixel size is in um
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

tic
disp('Inferring particle motion model...')
HMMParticles = track02TrainGHMM(RawParticles, FrameInfo, retrack, displayFigures);
toc

tic
disp('Simulating particle tracks...')
SimParticles = track03PredictParticlePaths(HMMParticles, FrameInfo, retrack, displayFigures);
toc

tic
disp('Stitching particle tracks...')
Particles = track04StitchTracks(SimParticles, FrameInfo, experimentType,...
                                useHistone, retrack, displayFigures);
toc 

disp('Adding QC fields...')
% Iterate over all channels and generate additional QC flags
maxCost = 3;
ncDistPrctile = 99.5;
for Channel = 1:nCh
  allDistanceVec = [Particles{Channel}.NucleusDist];  
  threshDist = prctile(allDistanceVec,ncDistPrctile);
  for p = 1:length(Particles{Channel})
    % flag particles anomalously far from their respective nuclei
    ncDistVec = Particles{Channel}(p).NucleusDist;

    % flag cases when particle is far away from nearest nucleus
    Particles{Channel}(p).ncDistFlags = ncDistVec>threshDist;
    
    % flag unlikely linkages
    Particles{Channel}(p).linkFlags = Particles{Channel}(p).linkCosts>maxCost;    
  end    
end


% for currentChannel = 1:nCh
%     
%     if ~isfield(Particles{currentChannel}, 'FrameApproved')
%         
%         for particle = 1:length(Particles{currentChannel})
%             Particles{currentChannel}(particle).FrameApproved = true(size(Particles{currentChannel}(particle).Frame));
%         end
%         
%     else
%         
%         for particle = 1:length(Particles{currentChannel})
%             
%             if isempty(Particles{currentChannel}(particle).FrameApproved)
%                 Particles{currentChannel}(particle).FrameApproved = true(size(Particles{currentChannel}(particle).Frame));
%             end
%             
%         end
%         
%     end
%     
%     
%     Particles = addPositionsToParticles(Particles, Spots, currentChannel);
%     
% end



% If we only have one channel, then convert SpotFilter and Particles to a standard structure.
if nCh == 1
    Particles = Particles{1};
end

end