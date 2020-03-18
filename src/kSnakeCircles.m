function [mask, ellipseFrame] = kSnakeCircles(image, PixelSize_um)

%parameters i've found to be broadly applicable
sigmaK= 8; 
sigmaSnake = sigmaK;
b = -.4; 
s = .1;

min_rad_um = 2; % set min and max acceptable area for nucleus segmentation
max_rad_um = 6; %this needs to be 6um for nc12. 4um for nc14
minArea_px = round(pi*(min_rad_um ./ PixelSize_um).^2);
maxArea_px = round(pi*(max_rad_um ./ PixelSize_um).^2);

areaFilter = [minArea_px, maxArea_px];

%preprocessing denoise
image = wiener2(image);


%assuming the right class is 3 here. sometimes that might be wrong.
%it's hard to decide which one is right automatically.
kmask = imsegkmeans(single(imgaussfilt(image,sigmaK)),3)==3;
%sometimes snakes destroys blobs. if it does, it'd be nice to add back in
%regions from kmask.

snakesFun = @(b, s, sigma) activecontour(imgaussfilt(image, sigma), kmask, 'Chan-Vese', 'ContractionBias', b, 'SmoothFactor', s);
kMaskRefined = snakesFun(b, s, sigmaK);
%same issue here. sometimes the watershed is inverted. hard to pick
%whether it should be or not automatically

mask = bwareafilt(wshed(kMaskRefined), areaFilter);

%fit with circles instead of convex hulls
[mask, ellipseFrame] = fitCirclesToNuclei(mask, kmask);


end