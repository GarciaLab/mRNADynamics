function schnitzcells = TrackNuclei_computervision(Prefix)
%AR 4/2020
load('ReferenceHist.mat', 'ReferenceHist');
thisExperiment = liveExperiment(Prefix);
pixelSize_um = thisExperiment.pixelSize_um;
nFrames = thisExperiment.nFrames;
hisVideoFile = [thisExperiment.preFolder, filesep, 'hisVideo.avi'];

if ~exist(hisVideoFile, 'file')
    hisMat = getHisMat(thisExperiment);
    exportTifStackToAvi(hisMat, hisVideoFile)
end


% memoMasker = memoize (@(x, y) kSnakeCircles(x, y));
% memoMasker.CacheSize = nFrames;

obj = setupSystemObjects(hisVideoFile);
tracks = initializeTracks(); % Create an empty array of tracks.
nextId = 1; % ID of the next track
% Detect moving objects, and track them across video frames.
frameIndex = 0;

for f = 1:nFrames
    
    frameIndex = frameIndex + 1;
    
    frame = readFrame(obj.reader);
    
    [centroids, radii, bboxes, mask] =...
        detectObjects(obj, frame, pixelSize_um);
    
    tracks = predictNewLocationsOfTracks(tracks);
    
    [assignments, unassignedTracks,...
        unassignedDetections] = ...
        detectionToTrackAssignment(tracks, centroids);
    
    tracks = updateAssignedTracks(tracks, assignments, bboxes,...
        centroids);
    
    tracks = updateUnassignedTracks(tracks, unassignedTracks);
    tracks = deleteLostTracks(tracks);
    
    tracks = createNewTracks(tracks, centroids, bboxes, radii);
    
    schnitzcells = displayTrackingResults(tracks, obj, frame,...
        mask, ReferenceHist);
    
end

end