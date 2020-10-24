    function particleTracks = updateAssignedParticleTracks(...
                              particleTracks, assignments, measurements, CurrentFrame,...
                              zOrig)
        
        numAssignedTracks = 0;
        if size(assignments, 2)==2
          numAssignedTracks = size(assignments, 1);
        end
        
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
            particleTracks(trackIdx).Index(end + 1) = detectionIdx;
            particleTracks(trackIdx).MeasurementVec(end + 1,:) = measurement;
            particleTracks(trackIdx).zPos(end + 1) = zOrig(detectionIdx);
            
%             if size(particleTracks(trackIdx).MeasurementVec,1)>1
%               if abs(particleTracks(trackIdx).MeasurementVec(end-1,1)-particleTracks(trackIdx).MeasurementVec(end,1)) > 20 && ...
%                   particleTracks(trackIdx).Frame(end-1)-particleTracks(trackIdx).Frame(end)==1
%                   error('wtf')
%               end
%             end
            % Update visibility.
            particleTracks(trackIdx).totalVisibleCount = ...
                particleTracks(trackIdx).totalVisibleCount + 1;
            particleTracks(trackIdx).consecutiveInvisibleCount = 0;
            
            % update last frame info
            particleTracks(trackIdx).lastFrame = max(particleTracks(trackIdx).Frame);  
            particleTracks(trackIdx).firstFrame = min(particleTracks(trackIdx).Frame);  
            
        end
        
    end
