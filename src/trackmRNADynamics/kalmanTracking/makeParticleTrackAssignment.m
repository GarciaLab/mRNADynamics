function [assignments, unassignedTracks, unassignedDetections] = ...
                  makeParticleTrackAssignment(particleTracks, measurements,...
                  maxCost,NewSpotNuclei,maxUnobservedFrames)

  nTracks = length(particleTracks);
  nDetections = size(measurements, 1);

  % Compute the cost of assigning each detection to each track.

  costArray = zeros(nTracks, nDetections);
  
  % NL: This returns the negative of  the log likelihood, with 
  % trajectories modeled as multidimensional Gaussian process
  extantNuclei = [particleTracks.Nucleus];
  cappedParticles = [particleTracks.consecutiveInvisibleCount]>=maxUnobservedFrames;
  for k = 1:nTracks
      costArray(k, :) = distance(particleTracks(k).kalmanFilter, measurements);
      costArray(k, NewSpotNuclei~=extantNuclei(k)) = Inf; % forbid matches between spots from different nuclei
  end
  costArray(cappedParticles,:) = Inf;
  
  % perform pairwise matching
  [assignments, unassignedTracks, unassignedDetections] = ...
                        assignDetectionsToTracks(costArray, 0.5*maxCost);

end