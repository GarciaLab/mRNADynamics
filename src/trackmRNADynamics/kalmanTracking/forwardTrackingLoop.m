function [forwardTracks, trackingInfo] = forwardTrackingLoop(forwardTracks, trackingInfo, ...
            kalmanOptions, Spots, Channel, CurrentFrame, schnitzcells, calibrateCost)
          
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
            schnitzcells,CurrentFrame,trackingInfo.useHistone);

    if kalmanOptions.useNuclei
      SpotMeasurements = [SpotMeasurements NucleusDistances];
    end
    
    % initialize to proper dimensions if empty
    if isempty(SpotMeasurements) && kalmanOptions.useNuclei
      SpotMeasurements = zeros(1,5);
    elseif isempty(SpotMeasurements)
      SpotMeasurements = zeros(1,4);
    end
    % Check if we're at the start of a new nuclear cycle
    continuedNCFlag = true;
    if CurrentFrame > 1
      continuedNCFlag = trackingInfo.ncVec(CurrentFrame-1)==trackingInfo.ncVec(CurrentFrame);
    end                        

    % predict positions of extant particles
    forwardTracks = predictParticleTrackLocations(forwardTracks);

    earlyFlags = false(size(forwardTracks));

    if false % NL: this chunk is meant to guess an optimal linking cost...not working well yet
        
        [assignments, unassignedTracks, unassignedDetections, costArray] = ...
                    makeParticleTrackAssignment(forwardTracks, SpotMeasurements(:,trackingInfo.trackingIndices), ...
                    100, NewSpotNuclei,...
                    [], [], earlyFlags);
                  
        % get indices of assignments
        linearIndices = sub2ind(size(costArray),assignments(:,1),assignments(:,2));
        matchCosts = costArray(linearIndices);
        pArray = NaN(100,100);
        for p = 1:size(pArray,2)
          for n = 1:size(pArray,1)
            matchCostsBoot = randsample(matchCosts,length(matchCosts),true);
            pArray(n,p) = prctile(matchCostsBoot,p);
          end
        end
        trackingInfo.matchCostMax(Channel) = mean(costVec(iter-2:iter-1));
        
    else
        % Perform cost-based matching
        [assignments, unassignedTracks, unassignedDetections] = ...
                    makeParticleTrackAssignment(forwardTracks, SpotMeasurements(:,trackingInfo.trackingIndices), ...
                    continuedNCFlag*trackingInfo.matchCostMaxForward(Channel), NewSpotNuclei,...
                    [], [], earlyFlags);
    end

    % make new entries for spots that were not assigned to existing
    % particles
    forwardTracks = makeNewTracks(forwardTracks, SpotMeasurements,...
                                     unassignedDetections, kalmanOptions,...
                                     CurrentFrame, NewSpotsZ, NewSpotNuclei, trackingInfo);

    % update existing tracks that had no match this frame
    forwardTracks = updateUnassignedParticleTracks(forwardTracks, unassignedTracks,...
              trackingInfo.maxUnobservedFrames(Channel));          

    % update tracks that matched with a new particle
    forwardTracks = updateAssignedParticleTracks(...
                            forwardTracks, assignments, SpotMeasurements, CurrentFrame,...
                            NewSpotsZ, trackingInfo);     

    % lastly, identify and "Cap" cases where there are too many tracks
    % per nucleus
    forwardTracks = cleanUpTracks(forwardTracks, trackingInfo.spotsPerNucleus(Channel), ...
      ~continuedNCFlag||CurrentFrame==trackingInfo.nFrames); 