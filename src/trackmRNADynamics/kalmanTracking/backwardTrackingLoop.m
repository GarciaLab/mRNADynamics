function [backwardTracks, trackingInfo] = backwardTrackingLoop(backwardTracks, trackingInfo, ...
            kalmanOptions, Spots, Channel, CurrentFrame, schnitzcells)

    % Get the positions of ALL spots (approved and disapproved)
    [NewSpotsX, NewSpotsY, NewSpotsZ, NewSpotsFluo] = SpotsXYZ(Spots{Channel}(CurrentFrame));

    % adjust Z position variable for stage movements
    NewSpotsZAdjusted = NewSpotsZ - trackingInfo.zPosStage(CurrentFrame);
    SpotMeasurements = [NewSpotsX', NewSpotsY', NewSpotsZAdjusted'];

    % rescale Fluorescence values
    NewSpotsFluoRescaled = NewSpotsFluo / kalmanOptions.fluoFactor;
    SpotMeasurements = [SpotMeasurements NewSpotsFluoRescaled'];

    % if we have nucleus info, assign each spot to a nucleus
    [NewSpotNuclei, NucleusDistances] = getNuclearAssigments(NewSpotsX,NewSpotsY,...
            schnitzcells, CurrentFrame, trackingInfo.useHistone);

    if false % NL: leave out nucleus info for now...
        SpotMeasurements = [SpotMeasurements NucleusDistances];
    end                

    % initialize to proper dimensions if empty
    if isempty(SpotMeasurements) && kalmanOptions.useNuclei
      SpotMeasurements = zeros(1,5);
    elseif isempty(SpotMeasurements)
      SpotMeasurements = zeros(1,4);
    end
    
    % get info about forward tracks taht were active during this frame
    activeArrayIndices = find(trackingInfo.activeArray(CurrentFrame,:)==1);
    activeSpotIndices = trackingInfo.indexArray(CurrentFrame,activeArrayIndices);
    activeParticleIndices = trackingInfo.assignmentArray(CurrentFrame,activeArrayIndices);
    assignedParticles = trackingInfo.assignmentArray(CurrentFrame,:);

    % Check if we're at the start of a new nuclear cycle
    continuedNCFlag = true;
    if CurrentFrame < trackingInfo.nFrames
      continuedNCFlag = trackingInfo.ncVec(CurrentFrame+1)==trackingInfo.ncVec(CurrentFrame);
    end      

    % predict positions of extant particles
    backwardTracks = predictParticleTrackLocations(backwardTracks);

    % get list of tracks that have upcoming assigned particles
    ffVec = NaN(size(backwardTracks));
    for p = 1:length(ffVec)
        ffVec(p) = min(trackingInfo.firstFrameVec(assignedParticles==p));
    end
    
    earlyFlags = ffVec<CurrentFrame;
    % Perform cost-based matching
    [assignments, unassignedTracks, unassignedDetections] = ...
                makeParticleTrackAssignment(backwardTracks, SpotMeasurements(:,trackingInfo.trackingIndices), ...
                continuedNCFlag*trackingInfo.matchCostMaxBackward(Channel), NewSpotNuclei,...
                activeSpotIndices, activeParticleIndices,earlyFlags);

    % update mapping vec     
    if ~isempty(assignments)
      for i = 1:size(assignments,1)          
          colIndex = trackingInfo.indexArray(CurrentFrame,:)==assignments(i,2);
          if all(isnan(trackingInfo.assignmentArray(:,colIndex)))
              trackingInfo.assignmentArray(:,colIndex) = assignments(i,1);
          end
      end     
    end
    newIndVec = length(backwardTracks)+1:length(backwardTracks)+length(unassignedDetections);
    for i = 1:length(unassignedDetections)    
        colIndex = trackingInfo.indexArray(CurrentFrame,:)==unassignedDetections(i);         
        trackingInfo.assignmentArray(:,colIndex) = newIndVec(i);
    end

    % make new entries for spots that were not assigned to existing
    % particles      
    backwardTracks = makeNewTracks(backwardTracks, SpotMeasurements,...
                                     unassignedDetections, kalmanOptions,...
                                     CurrentFrame, NewSpotsZ, NewSpotNuclei, trackingInfo);

    % update existing tracks that had no match this frame
    backwardTracks = updateUnassignedParticleTracks(backwardTracks, unassignedTracks, ...
            trackingInfo.maxUnobservedFrames(Channel));          

    % update tracks that matched with a new particle
    backwardTracks = updateAssignedParticleTracks(...
                            backwardTracks, assignments, SpotMeasurements,...
                            CurrentFrame, NewSpotsZ, trackingInfo);   


   % lastly, identify and "Cap" cases where there are too many tracks
    % per nucleus
    backwardTracks = cleanUpTracks(backwardTracks, ...
      trackingInfo.spotsPerNucleus(Channel), ~continuedNCFlag||CurrentFrame==1);  

