    function [tracks, nextId] = createNewTracks(tracks, measurements,...
        bboxes, unassignedDetections, nextId)
       
        %measurements. should be cated
        measurements = measurements(...
        unassignedDetections, :);

        bboxes = bboxes(unassignedDetections, :);
        
        for i = 1:size(measurements, 1)
            
            bbox = bboxes(i, :);
            measurement = measurements(i, :);
            
            % Create a Kalman filter object.
            
            InitialEstimateError = [50, 50, 50];
            MotionNoise = [25, 25, 25];
            MeasurementNoise = 50; %pixel variance of 100
            
            %     kalmanFilter = configureKalmanFilter(MotionModel,InitialLocation,InitialEstimateError,...
            %         MotionNoise,MeasurementNoise)
            kalmanFilter = configureKalmanFilter('ConstantAcceleration', ...
                measurement, InitialEstimateError , MotionNoise,  MeasurementNoise );
            
            % Create a new track.
            newTrack = struct(...
                'id', nextId, ...
                'bbox', bbox, ...
                'kalmanFilter', kalmanFilter, ...
                'age', 1, ...
                'totalVisibleCount', 1, ...
                'consecutiveInvisibleCount', 0);
            
            % Add it to the array of tracks.
            tracks(end + 1) = newTrack;
            
            % Increment the next id.
            nextId = nextId + 1;
        end
        
    end