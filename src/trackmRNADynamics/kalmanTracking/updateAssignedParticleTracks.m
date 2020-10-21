    function particleTracks = updateAssignedParticleTracks(...
                              particleTracks, assignments, measurements, CurrentFrame,...
                              zOrig, NewSpotNuclei)
        
        numAssignedTracks = size(assignments, 1);
        
        for i = 1:numAssignedTracks
            
            trackIdx = assignments(i, 1);
            detectionIdx = assignments(i, 2);
            measurement = measurements(detectionIdx, :);            
            
            % Correct the estimate of the object's location
            % using the new detection. This modifies
            % the current kalmanFilter object
            % in memory. 
            correct(particleTracks(trackIdx).kalmanFilter, measurement);          
            
            % Update track's age.
            particleTracks(trackIdx).age = particleTracks(trackIdx).age + 1;  
            
            % Update core particle parameters
            particleTracks(trackIdx).Frame(end + 1) = CurrentFrame;
            particleTracks(trackIdx).MeasurementVec(end + 1,:) = measurement;
%             particleTracks(trackIdx).Nucleus(end + 1) = NewSpotNuclei(detectionIdx);
            particleTracks(trackIdx).zPos(end + 1) = zOrig(detectionIdx);
            % Update visibility.
            particleTracks(trackIdx).totalVisibleCount = ...
                particleTracks(trackIdx).totalVisibleCount + 1;
            particleTracks(trackIdx).consecutiveInvisibleCount = 0;
            
        end
        
    end
