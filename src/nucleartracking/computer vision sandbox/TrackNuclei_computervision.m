function schnitzcells = TrackNuclei_computervision(Prefix)

load('ReferenceHist.mat', 'ReferenceHist');
thisExperiment = liveExperiment(Prefix);


hisVideoFile = [thisExperiment.preFolder, filesep, 'hisVideo.avi'];

if ~exist(hisVideoFile, 'file')
    hisMat = getHisMat(thisExperiment);    
    exportTifStackToAvi(hisMat, hisVideoFile)
end

schnitzcells = struct('cenx', [], 'ceny', [],...
    'frames', []);

pixelSize_um = thisExperiment.pixelSize_um;
nFrames = thisExperiment.nFrames;
% memoMasker = memoize (@(x, y) kSnakeCircles(x, y));
% memoMasker.CacheSize = nFrames;

% Ellipses = getEllipses(thisExperiment)
obj = setupSystemObjects(hisVideoFile);
tracks = initializeTracks(); % Create an empty array of tracks.
nextId = 1; % ID of the next track
% Detect moving objects, and track them across video frames.
frameIndex = 0;

% while hasFrame(obj.reader)
for f = 1:nFrames
    frameIndex = frameIndex + 1;
    frame = readFrame(obj.reader);
    [centroids, bboxes, mask] = detectObjects(frame);
    predictNewLocationsOfTracks();
    [assignments, unassignedTracks, unassignedDetections] = ...
        detectionToTrackAssignment();
    
    updateAssignedTracks();
    updateUnassignedTracks();
    deleteLostTracks();
    createNewTracks();
   
    displayTrackingResults();
    
end


function tracks = initializeTracks()
% create an empty array of tracks
tracks = struct(...
    'id', {}, ...
    'bbox', {}, ...
    'kalmanFilter', {}, ...
    'age', {}, ...
    'totalVisibleCount', {}, ...
    'consecutiveInvisibleCount', {});
end

function [centroids, bboxes, mask] = detectObjects(frame)
% Detect foreground.
% mask = obj.detector.step(frame);

[mask, ~] = kSnakeCircles(frame, pixelSize_um);
% [mask, ~] = memoMasker(frame, pixelSize_um);

mask = ~~mask;
% Perform blob analysis to find connected components.
[~, centroids, bboxes] = obj.blobAnalyser.step(mask);

%bbox is x y w h

end

function predictNewLocationsOfTracks()
for i = 1:length(tracks)
    bbox = tracks(i).bbox;
    
    % Predict the current location of the track.
    predictedCentroid = predict(tracks(i).kalmanFilter);
    
    % Shift the bounding box so that its center is at
    % the predicted location.
    predictedCentroid = int32(predictedCentroid) - bbox(3:4) / 2;
    tracks(i).bbox = [predictedCentroid, bbox(3:4)];
    
%     tracks(i).cenx = predictedCentroid(2);
%     tracks(i).ceny = predictedCentroid(1);
end
end

function [assignments, unassignedTracks, unassignedDetections] = ...
    detectionToTrackAssignment()

nTracks = length(tracks);
nDetections = size(centroids, 1);

% Compute the cost of assigning each detection to each track.
cost = zeros(nTracks, nDetections);
for i = 1:nTracks
    cost(i, :) = distance(tracks(i).kalmanFilter, centroids);
end

% Solve the assignment problem.
costOfNonAssignment = 20;
[assignments, unassignedTracks, unassignedDetections] = ...
    assignDetectionsToTracks(cost, costOfNonAssignment);
end

function updateAssignedTracks()
numAssignedTracks = size(assignments, 1);
for i = 1:numAssignedTracks
    trackIdx = assignments(i, 1);
    detectionIdx = assignments(i, 2);
    centroid = centroids(detectionIdx, :);
    bbox = bboxes(detectionIdx, :);
    
    % Correct the estimate of the object's location
    % using the new detection.
    correct(tracks(trackIdx).kalmanFilter, centroid);
    
    % Replace predicted bounding box with detected
    % bounding box.
    tracks(trackIdx).bbox = bbox;
    
    % Update track's age.
    tracks(trackIdx).age = tracks(trackIdx).age + 1;
    
    % Update visibility.
    tracks(trackIdx).totalVisibleCount = ...
        tracks(trackIdx).totalVisibleCount + 1;
    tracks(trackIdx).consecutiveInvisibleCount = 0;
end
end

function updateUnassignedTracks()
for i = 1:length(unassignedTracks)
    ind = unassignedTracks(i);
    tracks(ind).age = tracks(ind).age + 1;
    tracks(ind).consecutiveInvisibleCount = ...
        tracks(ind).consecutiveInvisibleCount + 1;
end
end

function deleteLostTracks()
if isempty(tracks)
    return;
end

invisibleForTooLong = 100;
ageThreshold = 1; %8

% Compute the fraction of the track's age for which it was visible.
ages = [tracks(:).age];
totalVisibleCounts = [tracks(:).totalVisibleCount];
visibility = totalVisibleCounts ./ ages;

% Find the indices of 'lost' tracks.
lostInds = (ages < ageThreshold & visibility < 0.6) | ...
    [tracks(:).consecutiveInvisibleCount] >= invisibleForTooLong;

% Delete lost tracks.
tracks = tracks(~lostInds);
end


function createNewTracks()
centroids = centroids(unassignedDetections, :);
bboxes = bboxes(unassignedDetections, :);

for i = 1:size(centroids, 1)
    
    centroid = centroids(i,:);
    bbox = bboxes(i, :);
    
    % Create a Kalman filter object.
    InitialEstimateError = [50, 50, 50];
    MotionNoise = [25, 25, 25];
    MeasurementNoise = 50; %pixel variance of 100
    
%     kalmanFilter = configureKalmanFilter(MotionModel,InitialLocation,InitialEstimateError,...
%         MotionNoise,MeasurementNoise)
    kalmanFilter = configureKalmanFilter('ConstantAcceleration', ...
        [centroid, size(bbox, 1)], InitialEstimateError , MotionNoise,  MeasurementNoise );
    A_n = [1 1 .5
                0 1 1
                0   0 1];
    %one block for x, one for y, one for radius       
%     MotionModel = blkdiag(A_n, A_n, A_n);
%     kalmanFilter = vision.kalmanFilter(MotionModel,[centroid, size(bbox, 1)],InitialEstimateError,...
%         MotionNoise,MeasurementNoise);
% 
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

function displayTrackingResults()
% Convert the frame and the mask to uint8 RGB.

mask = histeq(uint8(mask), ReferenceHist);

minVisibleCount = 0;

if ~isempty(tracks)
    
    % Noisy detections tend to result in short-lived tracks.
    % Only display tracks that have been visible for more than
    % a minimum number of frames.
    reliableTrackInds = ...
        [tracks(:).totalVisibleCount] > minVisibleCount;
    reliableTracks = tracks(reliableTrackInds);
    
    % Display the objects. If an object has not been detected
    % in this frame, display its predicted bounding box.
    if ~isempty(reliableTracks)
        % Get bounding boxes.
        bboxes = cat(1, reliableTracks.bbox);
        cenxs = (bboxes(:, 1) + bboxes(:, 3))/2;
        cenys = (bboxes(:, 2) + bboxes(:, 4))/2;
        % Get ids.
        ids = int32([reliableTracks(:).id]);
        
        for i = 1:length(ids)
            if length(schnitzcells) < i
                
            schnitzcells(i).cenx(1) = cenxs(i);
            schnitzcells(i).ceny(1) = cenys(i);
            schnitzcells(i).frames(1) = frameIndex;
            
            else
                
            schnitzcells(i).cenx(end+1) = cenxs(i);
            schnitzcells(i).ceny(end+1) = cenys(i);
            schnitzcells(i).frames(end+1) = frameIndex;
            
            end
        end
        
        % Create labels for objects indicating the ones for
        % which we display the predicted rather than the actual
        % location.
        labels = cellstr(int2str(ids'));
        predictedTrackInds = ...
            [reliableTracks(:).consecutiveInvisibleCount] > 0;
        isPredicted = cell(size(labels));
        isPredicted(predictedTrackInds) = {' predicted'};
        labels = strcat(labels, isPredicted);
        
        % Draw the objects on the frame.
        frame = insertObjectAnnotation(frame, 'rectangle', ...
            bboxes, labels);
        frame = rgb2gray(frame);
        % Draw the objects on the mask.
        mask = insertObjectAnnotation(mask, 'rectangle', ...
            bboxes, labels);
        mask = rgb2gray(mask);
    end
end

% Display the mask and the frame.
obj.maskPlayer.step(mask);
obj.videoPlayer.step(frame);

end

end

function obj = setupSystemObjects(hisVideoFile)
% 
obj.reader = VideoReader(hisVideoFile);
obj.maskPlayer = vision.VideoPlayer('Position', [255 205 415 300]);
obj.videoPlayer = vision.VideoPlayer('Position', [726 205 415 300]);

obj.detector = vision.ForegroundDetector('NumGaussians', 8, ...
            'NumTrainingFrames', 40, 'MinimumBackgroundRatio', 0.4);

obj.blobAnalyser = vision.BlobAnalysis('BoundingBoxOutputPort', true, ...
            'AreaOutputPort', true, 'CentroidOutputPort', true, ...
            'MinimumBlobArea', 400); 
        
end