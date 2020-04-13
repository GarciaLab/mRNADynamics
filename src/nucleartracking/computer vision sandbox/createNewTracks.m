    function [tracks, nextId] = createNewTracks(tracks, measurements,...
        bboxes, unassignedDetections, nextId)
       
        %ceny, cenx, smaj, smin, angle
        nStates = 5;
    
        measurements = measurements(...
        unassignedDetections, :);

        bboxes = bboxes(unassignedDetections, :);
        
        for i = 1:size(measurements, 1)
            
            bbox = bboxes(i, :);
            measurement = measurements(i, :);
            
            % Create a Kalman filter object.
            centroidError_pxSq = 4^2;
            radiusError_pxSq = 5^2;
            angleError_radSq = (pi/2)^2;
            
            InitialEstimateError = [1 1 1]*1e5;
            MotionNoise = [centroidError_pxSq, centroidError_pxSq ,...
                radiusError_pxSq];
            MeasurementNoise = radiusError_pxSq;
            
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
                'idxHistory', nextId,...
                'totalVisibleCount', 1, ...
                'consecutiveInvisibleCount', 0);
            
            % Add it to the array of tracks.
            tracks(end + 1) = newTrack;
            
            % Increment the next id.
            nextId = nextId + 1;
        end
        
    end