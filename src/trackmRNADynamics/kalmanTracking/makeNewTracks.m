function particleTracks = makeNewTracks(particleTracks, measurements,...
                                       unassignedDetections, kalmanOptions,...
                                       CurrentFrame, zOrig, NewSpotNuclei, trackingInfo)
           
    % get specs for unassigned particles
    measurements = measurements(unassignedDetections, :);

    for i = 1:size(measurements, 1)

        measurement = measurements(i, :);

        % Create a Kalman filter object            

        kalmanFilter = configureKalmanFilter(kalmanOptions.type, ...
            measurement(trackingInfo.trackingIndices), kalmanOptions.InitialError, kalmanOptions.MotionNoise,  kalmanOptions.MeasurementNoise);

        % Create a new track.
        newTrack = struct(...
            'Index', unassignedDetections(i), ...  
            'Frame', CurrentFrame, ... 
            'MeasurementVec', measurement, ...                 
            'zPos', zOrig(unassignedDetections(i)), ... % NL: need to keep track of unadjusted Z, even though we don't use for tracking 
            'Nucleus', NewSpotNuclei(unassignedDetections(i)), ...
            'kalmanFilter', kalmanFilter, ...
            'age', 1, ...                
            'totalVisibleCount', 1, ...
            'consecutiveInvisibleCount', 0,...
            'firstFrame', CurrentFrame, ...
            'lastFrame', CurrentFrame, ...
            'cappedFlag', false, ...
            'duplicateFlag',false);

        % Add it to the array of tracks.
        particleTracks(end + 1) = newTrack;

    end
        
end