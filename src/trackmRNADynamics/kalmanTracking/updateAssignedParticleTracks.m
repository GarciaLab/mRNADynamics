    function tracks= updateAssignedParticleTracks(...
                            tracks, assignments, bboxes, measurements,CurrentFrame)
        
        numAssignedTracks = size(assignments, 1);
        
        for i = 1:numAssignedTracks
            
            trackIdx = assignments(i, 1);
            detectionIdx = assignments(i, 2);
            measurement = measurements(detectionIdx, :);
            bbox = bboxes(detectionIdx, :);
            
            % Correct the estimate of the object's location
            % using the new detection. This modifies
            %the current kalmanFilter object
            %in memory. 
            correct(tracks(trackIdx).kalmanFilter, measurement);

            % Replace predicted bounding box with detected
            % bounding box.
            tracks(trackIdx).bbox = bbox;
            
            % Update track's age.
            tracks(trackIdx).age = tracks(trackIdx).age + 1;            
            % Update track's history
            tracks(trackIdx).idxHistory=...
                [ tracks(trackIdx).idxHistory, trackIdx];
           
            % Update visibility.
            tracks(trackIdx).totalVisibleCount = ...
                tracks(trackIdx).totalVisibleCount + 1;
            tracks(trackIdx).consecutiveInvisibleCount = 0;
            
        end
        
    end
