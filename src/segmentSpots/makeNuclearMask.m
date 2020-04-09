function nuclearMask = makeNuclearMask(ellipseFrame, imageDims, radiusScale)

%make a mask from the ellipses structure. this is currently used for
%segmenting loci in segmentSpots.

if nargin < 3
    radiusScale = 1.3; %be more lenient with the size of ellipses in the nuclear mask
                                %so spots don't get excluded inappropriately
end

yDim = imageDims(1);
xDim = imageDims(2);

nuclearMask = false(yDim, xDim);

nEllipses = size(ellipseFrame, 1);

for n = 1:nEllipses
    h = images.roi.Ellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],...
        'SemiAxes', radiusScale*[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
        'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)),'StripeColor','m');
    nuclearMask = nuclearMask + poly2mask(h.Vertices(:, 1), h.Vertices(:, 2), yDim, xDim);
end

