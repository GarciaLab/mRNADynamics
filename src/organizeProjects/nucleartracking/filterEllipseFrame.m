function ellipseFrame = filterEllipseFrame(ellipseFrame, imageDims)

%validate sizes. the ellipse masker handles
%very large objects poorly, to say the least.

yDim = imageDims(1);
xDim = imageDims(2);

largeAxisIndex = ellipseFrame(:, 3) > max([xDim, yDim])...
    | ellipseFrame(3) > max([xDim, yDim]);
for k = 1:length(largeAxisIndex)
    if largeAxisIndex(k)
        ellipseFrame(k, :) = zeros(1, size(ellipseFrame, 2)); 
    end
end