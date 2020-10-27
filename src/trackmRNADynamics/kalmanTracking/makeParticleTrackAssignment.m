function [assignments, unassignedTracks, unassignedDetections, costArray] = ...
                  makeParticleTrackAssignment(particleTracks, measurements,...
                  maxCost,NewSpotNuclei, activeSpotIndices,...
                  activeParticleIndices, earlyFlags, trackingOptions) 
  
  % NL: This returns the negative of  the log likelihood, with 
  % trajectories modeled as multidimensional Gaussian process
  particleNuclei = [particleTracks.Nucleus];  
  cappedParticles = [particleTracks.cappedFlag];
  duplicateParticles = [particleTracks.duplicateFlag];     
  
  mappedParticles = ismember(1:length(particleTracks),activeParticleIndices); 
  eligibleParticles = ~cappedParticles & ~mappedParticles & ~earlyFlags & ~duplicateParticles;
  mappedDetectionIDs = activeSpotIndices(~isnan(activeParticleIndices));
  mappedParticleIDs = activeParticleIndices(~isnan(activeParticleIndices));

  % get list of unique nuclei for new detections
  nucleusIndex = unique(NewSpotNuclei);
  
  % assign pre-determined matches (this only applies for backwards
  % tracking)
  assignments = [mappedParticleIDs' mappedDetectionIDs'];
  
  % initialize list of unassigned tracks with capped particles  
  unassignedTracks = find(~cappedParticles)';  
  unassignedDetections = [];
  
  for n = 1:length(nucleusIndex)
    
      detectionIDs = find(NewSpotNuclei==nucleusIndex(n));
      detectionIDs = detectionIDs(~ismember(detectionIDs,mappedDetectionIDs));
      particleIDs = find(particleNuclei==nucleusIndex(n)&eligibleParticles);
      
      % Compute the cost of assigning each detection to each track.
      costArray = zeros(length(particleIDs), length(detectionIDs));

      for k = 1:length(particleIDs)
          costArray(k, :) = distance(particleTracks(particleIDs(k)).kalmanFilter, measurements(detectionIDs,:));          
      end      
            
      % perform pairwise matching
      [assignmentsTemp, unassignedTracksTemp, unassignedDetectionsTemp] = ...
                            assignDetectionsToTracks(costArray, 0.5*maxCost);

      % map back to full arrays
      assignments = [assignments ; [particleIDs(assignmentsTemp(:,1))',detectionIDs(assignmentsTemp(:,2))']];
      unassignedTracks = [unassignedTracks ; particleIDs(unassignedTracksTemp)'];
      unassignedDetections = [unassignedDetections ; detectionIDs(unassignedDetectionsTemp)'];
      
  end                           

end