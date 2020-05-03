function [schnitzcells, Ellipses] = TrackNuclei_computervision(Prefix, frameRange)

arguments
   Prefix char 
   frameRange (1, 2) double = [1 inf]
end

%AR 4/2020
load('ReferenceHist.mat', 'ReferenceHist');
liveExperiment = LiveExperiment(Prefix);
pixelSize_um = liveExperiment.pixelSize_um;

%if no frameRange upper limit was passed, change the upper limit 
%to the total number of frames in the movie
nFrames = liveExperiment.nFrames;
if isinf(frameRange(2))
    frameRange(2) = nFrames;
end


hisVideoFile = [liveExperiment.preFolder, filesep, 'hisVideo.avi'];
schnitzcellsFile = [liveExperiment.resultsFolder, filesep,Prefix,'_lin.mat'];
ellipsesFile = [liveExperiment.resultsFolder, filesep 'Ellipses.mat'];
schnitzcells = struct('cenx', [], 'ceny', [],...
    'frames', [], 'smaj', [], 'smin', [],...
    'orientationAngle', []);

if ~exist(hisVideoFile, 'file')
    hisMat = getHisMat(liveExperiment);
    exportTifStackToAvi(hisMat, hisVideoFile)
end

playerObj = setupSystemObjects(hisVideoFile);
tracks = initializeTracks(); % Create an empty array of tracks.
nextId = 1; % ID of the next track


% Detect moving objects, and track them across video frames.
for frameIndex = frameRange
        
    frameImage = readFrame(playerObj.reader);
    %% Measurement
    
    %Segment the nuclei to create measurements
    %that will be fed into Kalman filter.
    [measurements, bboxes, mask] =...
        detectObjects(frameImage, pixelSize_um, nFrames);
    
    %% Prediction
    
    tracks = predictNewLocationsOfTracks(tracks);
    
    %% Updating tracks
    
    [assignments, unassignedTracks,...
        unassignedDetections] = ...
        detectionToTrackAssignment(tracks, measurements);
    
    
    tracks = updateAssignedTracks(tracks, assignments, bboxes,...
        measurements);
    
    
    tracks = updateUnassignedTracks(tracks, unassignedTracks);
    
    tracks = deleteLostTracks(tracks);
    
    [tracks, nextId] = createNewTracks(tracks, measurements,...
        bboxes, unassignedDetections, nextId);
    
    schnitzcells = displayTrackingResults(tracks, playerObj, frameImage,...
        mask, ReferenceHist, schnitzcells, frameIndex);
    
end


[Ellipses, schnitzcells] = makeEllipsesFromSchnitzcells(schnitzcells, nFrames);

save2(ellipsesFile, Ellipses);
save2(schnitzcellsFile, schnitzcells);
% 
% schnitzcells = integrateSchnitzFluo(Prefix, schnitzcells,...
%     getFrameInfo(liveExperiment));
% 
% %%
% try
% expandedAnaphaseFrames = [zeros(1,8),liveExperiment.anaphaseFrames'];
% [schnitzcells, Ellipses] = breakUpSchnitzesAtMitoses(...
%     schnitzcells, Ellipses, expandedAnaphaseFrames, nFrames);
% save2(ellipsesFile, Ellipses);
% save2(schnitzcellsFile, schnitzcells);
% catch, warning('did not break up at mitoses');
% end
%%

%validate tracks
% figure;
% for k = 1:length(schnitzcells)
%     hold on
%     plot(schnitzcells(k).frames, schnitzcells(k).ceny, '-x')
%     plot(schnitzcells(k).frames, schnitzcells(k).len, '-x')
% end

end