function [backwardTracks, trackingOptions] = backwardTrackingLoop(backwardTracks, trackingOptions, ...
            kalmanOptions, Spots, Channel, CurrentFrame, schnitzcells)

    % Get the positions of ALL spots (approved and disapproved)
    [NewSpotsX, NewSpotsY, NewSpotsZ, NewSpotsFluo] = SpotsXYZ(Spots{Channel}(CurrentFrame),trackingOptions.use3DInfo);

    % adjust Z position variable for stage movements
    NewSpotsZAdjusted = NewSpotsZ - trackingOptions.zPosStage(CurrentFrame);
    SpotMeasurements = [NewSpotsX', NewSpotsY', NewSpotsZAdjusted'];

    % rescale Fluorescence values
    NewSpotsFluoRescaled = NewSpotsFluo / kalmanOptions.fluoFactor;
    SpotMeasurements = [SpotMeasurements NewSpotsFluoRescaled'];

    % if we have nucleus info, assign each spot to a nucleus
    [NewSpotNuclei, NucleusDistances] = getNuclearAssigments(NewSpotsX,NewSpotsY,...
            schnitzcells, CurrentFrame, trackingOptions.useHistone);

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
    activeArrayIndices = find(trackingOptions.activeArray(CurrentFrame,:)==1);
    activeSpotIndices = trackingOptions.indexArray(CurrentFrame,activeArrayIndices);
    activeParticleIndices = trackingOptions.assignmentArray(CurrentFrame,activeArrayIndices);
    assignedParticles = trackingOptions.assignmentArray(CurrentFrame,:);

    % Check if we're at the start of a new nuclear cycle
    continuedNCFlag = true;
    if CurrentFrame < trackingOptions.nFrames
      continuedNCFlag = trackingOptions.ncVec(CurrentFrame+1)==trackingOptions.ncVec(CurrentFrame);
    end      

    % predict positions of extant particles
    backwardTracks = predictParticleTrackLocations(backwardTracks);

    % get list of tracks that have upcoming assigned particles
    ffVec = NaN(size(backwardTracks));
    for p = 1:length(ffVec)
        ffVec(p) = min(trackingOptions.firstFrameVec(assignedParticles==p));
    end
    
    earlyFlags = ffVec<CurrentFrame;
    % Perform cost-based matching
    [assignments, unassignedTracks, unassignedDetections] = ...
                makeParticleTrackAssignment(backwardTracks, SpotMeasurements(:,trackingOptions.trackingIndices), ...
                continuedNCFlag*trackingOptions.matchCostMaxBackward(Channel), NewSpotNuclei,...
                activeSpotIndices, activeParticleIndices,earlyFlags);

    % update mapping vec     
    if ~isempty(assignments)
      for i = 1:size(assignments,1)          
          colIndex = trackingOptions.indexArray(CurrentFrame,:)==assignments(i,2);
          if all(isnan(trackingOptions.assignmentArray(:,colIndex)))
              trackingOptions.assignmentArray(:,colIndex) = assignments(i,1);
          end
      end     
    end
    newIndVec = length(backwardTracks)+1:length(backwardTracks)+length(unassignedDetections);
    for i = 1:length(unassignedDetections)    
        colIndex = trackingOptions.indexArray(CurrentFrame,:)==unassignedDetections(i);         
        trackingOptions.assignmentArray(:,colIndex) = newIndVec(i);
    end

    % make new entries for spots that were not assigned to existing
    % particles      
    backwardTracks = makeNewTracks(backwardTracks, SpotMeasurements,...
                                     unassignedDetections, kalmanOptions,...
                                     CurrentFrame, NewSpotsZ, NewSpotNuclei, trackingOptions);

    % update existing tracks that had no match this frame
    backwardTracks = updateUnassignedParticleTracks(backwardTracks, unassignedTracks, ...
            trackingOptions.maxUnobservedFrames(Channel));          

    % update tracks that matched with a new particle
    backwardTracks = updateAssignedParticleTracks(...
                            backwardTracks, assignments, SpotMeasurements,...
                            CurrentFrame, NewSpotsZ, trackingOptions);   


   % lastly, identify and "Cap" cases where there are too many tracks
    % per nucleus
    backwardTracks = cleanUpTracks(backwardTracks, ...
      trackingOptions.spotsPerNucleus(Channel), ~continuedNCFlag||CurrentFrame==1);  

