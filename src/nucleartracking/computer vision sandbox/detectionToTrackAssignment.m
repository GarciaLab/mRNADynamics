    function [assignments, unassignedTracks, unassignedDetections] = ...
            detectionToTrackAssignment(tracks, measurements)
        
        nTracks = length(tracks);
        nDetections = size(measurements, 1); 
        
        % Compute the cost of assigning each detection to each track.
        cost = zeros(nTracks, nDetections);
        for i = 1:nTracks
            cost(i, :) = distance(tracks(i).kalmanFilter, measurements);
        end
        
        % Solve the assignment problem.
        costOfNonAssignment = 20;
        [assignments, unassignedTracks, unassignedDetections] = ...
            assignDetectionsToTracks(cost, costOfNonAssignment);
    end