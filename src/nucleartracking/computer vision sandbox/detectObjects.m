    function [centroids,  radii, bboxes, mask] = detectObjects(obj, frame, pixelSize_um)
        % Detect foreground.
        % mask = obj.detector.step(frame);
        
        [mask, ellipseFrame] = kSnakeCircles(frame, pixelSize_um);
        % [mask, ~] = memoMasker(frame, pixelSize_um);
        mask = ~~mask; %binarize the mask
        centroids = ellipseFrame(:, 1:2);
        radii = (1/2) * mean(ellipseFrame(:, 3:4), 2);
        
        % Perform blob analysis to find connected components.
        [~, ~, bboxes] = obj.blobAnalyser.step(mask);
        %bbox is x y w h
        
    end

