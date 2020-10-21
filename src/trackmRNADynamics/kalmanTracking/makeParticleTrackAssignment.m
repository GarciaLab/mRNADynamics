function [assignments, unassignedTracks, unassignedDetections] = ...
    makeParticleTrackAssignment(tracks, measurements)

  nTracks = length(tracks);
  nDetections = size(measurements, 1);

  % Compute the cost of assigning each detection to each track.

  cost = zeros(nTracks, nDetections);
  % NL: This returns the negative of  the log likelihood, with 
  % trajectories modeled as multidimensional Gaussian process
  for k = 1:nTracks
      cost(k, :) = distance(tracks(k).kalmanFilter, measurements);
  end

  %AR- i don't know how to properly adjust
  %this number.
  costOfNonAssignment = -log(1e-3);%nanmin([prctile(min(cost),95) 20]);
  [assignments, unassignedTracks, unassignedDetections] = ...
      assignDetectionsToTracks(cost, costOfNonAssignment);

end