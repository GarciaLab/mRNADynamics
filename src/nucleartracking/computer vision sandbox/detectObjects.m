function [measurements, bboxes, mask] =...
    detectObjects(frame, pixelSize_um, nFrames)

%memoized for quicker debugging
%so we won't have to calculate this on
%repeated tracking calls
%         kSnakeCircles = memoize(@kSnakeCircles);
%         kSnakeCircles.CacheSize = nFrames;
[mask, ellipseFrame] = kSnakeCircles(frame, pixelSize_um);

mask = ~~mask; %binarize the mask
%         radii = (1/2) * mean(ellipseFrame(:, 3:4), 2);

if any(mask(:))
    %ceny, cenx, smaj, smin, angle
    measurements = ellipseFrame(:, 1:5);
else
    measurements = [];
end


boundingBoxCell = {};
for n = 1:size(ellipseFrame, 1)
    h = images.roi.Ellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
        'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)),'StripeColor','m');
    bw = createMask(h, size(frame, 1), size(frame, 2));
    props = regionprops(bw, 'BoundingBox');
    boundingBoxCell{n} = props.BoundingBox;
end

assert(size(measurements, 1) == numel(boundingBoxCell));

bboxes = [];
for k = 1:numel(boundingBoxCell)
    bboxes(k, :) = boundingBoxCell{k};
end
bboxes = int32(round(bboxes));

end

