    function tracks = createNewTracks(tracks, centroids,...
        bboxes, radii, unassignedDetections, nextId)
       
        centroids = centroids(unassignedDetections, :);
        bboxes = bboxes(unassignedDetections, :);
        
        for i = 1:size(centroids, 1)
            
            centroid = centroids(i,:);
            radius = radii(i);
            bbox = bboxes(i, :);
            measurement = [centroid, radius]; 
            
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