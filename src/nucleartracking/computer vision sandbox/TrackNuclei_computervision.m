function schnitzcells = TrackNuclei_computervision(Prefix)
%AR 4/2020
load('ReferenceHist.mat', 'ReferenceHist');
thisExperiment = liveExperiment(Prefix);
pixelSize_um = thisExperiment.pixelSize_um;
nFrames = thisExperiment.nFrames;
hisVideoFile = [thisExperiment.preFolder, filesep, 'hisVideo.avi'];
schnitzcells = struct('cenx', [], 'ceny', [],...
    'frames', [], 'len', []);

if ~exist(hisVideoFile, 'file')
    hisMat = getHisMat(thisExperiment);
    exportTifStackToAvi(hisMat, hisVideoFile)
end



obj = setupSystemObjects(hisVideoFile);
tracks = initializeTracks(); % Create an empty array of tracks.
nextId = 1; % ID of the next track
% Detect moving objects, and track them across video frames.
frameIndex = 0;

for f = 1:nFrames
    
    frameIndex = frameIndex + 1;
    
    frame = readFrame(obj.reader);
%% Measurement

    %Segment the nuclei to create measurements
    %that will be fed into Kalman filter. 
    [centroids, radii, bboxes, mask] =...
        detectObjects(frame, pixelSize_um, nFrames);
    measurements = [centroids, radii]; 
    
%% Prediction

    tracks = predictNewLocationsOfTracks(tracks);

%% Correction

    [assignments, unassignedTracks,...
        unassignedDetections] = ...
        detectionToTrackAssignment(tracks, measurements);
    
    tracks = updateAssignedTracks(tracks, assignments, bboxes,...
        measurements);
    
    tracks = updateUnassignedTracks(tracks, unassignedTracks);
    
    tracks = deleteLostTracks(tracks);
    
    [tracks, nextId] = createNewTracks(tracks, measurements,...
        bboxes, unassignedDetections, nextId);
    
    schnitzcells = displayTrackingResults(tracks, obj, frame,...
        mask, ReferenceHist, schnitzcells, frameIndex);
    
end

end