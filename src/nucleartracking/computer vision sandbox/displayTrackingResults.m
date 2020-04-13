function schnitzcells = displayTrackingResults(tracks, obj, frame,...
    mask, ReferenceHist, schnitzcells, frameIndex)

% Convert the frame and the mask to uint8 RGB.
% from logical.
mask = histeq(uint8(mask), ReferenceHist);

minVisibleCount = 0;

if ~isempty(tracks)
    
    % Noisy detections tend to result in short-lived tracks.
    % Only display tracks that have been visible for more than
    % a minimum number of frames.
%     reliableTrackInds = ...
%         [tracks(:).totalVisibleCount] > minVisibleCount;
%     reliableTracks = tracks(reliableTrackInds);
    reliableTracks = tracks;
    
    % Display the objects. If an object has not been detected
    % in this frame, display its predicted bounding box.
    if ~isempty(reliableTracks)
        % Get bounding boxes: [
        %[topleft_x, topleft_y,  width,  height]
        bboxes = cat(1, reliableTracks.bbox);
     
        % Get ids.
        ids = int32([reliableTracks(:).id]);
        
        schnitzcells = genSchnitzCellsFromIDs(schnitzcells, tracks,...
           frameIndex);
        
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