    function [centroids,  radii, bboxes, mask] =...
        detectObjects(frame, pixelSize_um, nFrames)
                
        %memoized for quicker debugging 
        %so we won't have to calculate this on
        %repeated tracking calls
        kSnakeCircles = memoize(@kSnakeCircles);
        kSnakeCircles.CacheSize = nFrames;
        [mask, ellipseFrame] = kSnakeCircles(frame, pixelSize_um);

        mask = ~~mask; %binarize the mask
        centroids = ellipseFrame(:, 1:2);
        radii = (1/2) * mean(ellipseFrame(:, 3:4), 2);
        
        %this code is silly, but i don't know how to 
        %shape the bbox matrix more elegantly.
        props = regionprops(mask, 'BoundingBox');
        a = {props.BoundingBox}; 
        for k = 1:numel(a)
            bboxes(k, :) = a{k};
        end
        bboxes = int32(round(bboxes));
                          
    end

