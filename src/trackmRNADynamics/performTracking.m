function Particles =...
    performTracking(Particles, schnitzcells, NCh, Spots, app,...
    PreProcPath, Prefix, UseHistone, ParticlesFig,...
    SpotsChannel, NucleiFig, particlesAxes, nucAxes,...
    Ellipses, PixelSize, SearchRadiusMicrons, ExperimentType,...
    FrameInfo, retrack, displayFigures)

% [Particles] = track02KalmanTesting(...
%     FrameInfo, Spots, NCh, PixelSize, SearchRadiusMicrons, retrack, displayFigures)
  
tic
disp('Performing intitial particle linking...')
RawParticles = track01ParticleProximity(...
    FrameInfo, Spots, schnitzcells, NCh, PixelSize, SearchRadiusMicrons, UseHistone, retrack, displayFigures);
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
Particles = track04StitchTracks(...
                          SimParticles, FrameInfo, ExperimentType, UseHistone, retrack, displayFigures);
toc                        
disp('Adding QC fields...')
% Iterate over all channels and generate additional QC flags
maxCost = 3;
ncDistPrctile = 99.5;
for Channel = 1:NCh
  allDistanceVec = [Particles{Channel}.NucleusDist];  
  threshDist = prctile(allDistanceVec,ncDistPrctile);
  for p = 1:length(Particles{Channel})
    % flag particles anomalously far from their respective nuclei
    ncDistVec = Particles{Channel}(p).NucleusDist;
%     frameVec = Particles{Channel}(p).Frame;
%     ncID = Particles{Channel}(p).NucleusID;
%     cellNo = schnitzcells(ncID).cellno;
%     meanRadiusVec = NaN(size(frameVec));
%     for f = 1:length(frameVec) 
%       meanRadiusVec(f) = mean(Ellipses{frameVec(f)}(cellNo(f),[3 4]));
%     end
    % flag cases when particle is far away from nearest nucleus
    Particles{Channel}(p).ncDistFlags = ncDistVec>threshDist;
    
    % flag unlikely linkages
    Particles{Channel}(p).linkFlags = Particles{Channel}(p).linkCosts>maxCost;    
  end    
end


% for currentChannel = 1:NCh
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
if NCh == 1
    Particles = Particles{1};
end

end