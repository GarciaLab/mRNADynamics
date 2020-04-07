function [mask, ellipseFrame] = kSnakeCircles(image,...
    PixelSize_um, varargin)

%parameters i've found to be broadly applicable
sigmaK_um = .85;
sigmaK_px = sigmaK_um / PixelSize_um;
mu = .1; %weight of length term for chen vese  algorithm. honestly don't know what this controls
min_rad_um = 2; % set min and max acceptable area for nucleus segmentation
max_rad_um = 6; %this needs to be 6um for nc12. 4um for nc14
nIterSnakes = 100;
% b = -.4; 
% s = .1;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for k = 1:2:(numel(varargin)-1)
    if k ~= numel(varargin)
        eval([varargin{k} '=varargin{k+1};']);
    end
end


minArea_px = round(pi*(min_rad_um ./ PixelSize_um).^2);
maxArea_px = round(pi*(max_rad_um ./ PixelSize_um).^2);

areaFilter = [minArea_px, maxArea_px];

%preprocessing denoise
image = wiener2(image);


%assuming the right class is 3 here. sometimes that might be wrong.
%it's hard to decide which one is right automatically.
%choose kMask does this. 
kLabel= imsegkmeans(single(imgaussfilt(image,sigmaK_px)),3);

kMask = kLabel == chooseKLabel(kLabel);

%sometimes snakes destroys blobs. if it does, it'd be nice to add back in
%regions from kMask.

% snakesFun = @(b, s, sigma) activecontour(imgaussfilt(image, sigma), kMask, 'Chan-Vese', 'ContractionBias', b, 'SmoothFactor', s);
% snakesFun = @(b, s, sigma) gather( chenvese( ...
%     imgaussfilt(gpuArray(image), sigma), 100, gpuArray(kMask))); %100 iterations
% kMaskRefined = snakesFun(b, s, sigmaK);


kMaskRefined= gather( chenvese( ...
    imgaussfilt( gpuArray(image), sigmaK_px),...
    gpuArray(kMask), nIterSnakes, mu, 'chan') );

%same issue here. sometimes the watershed is inverted. hard to pick
%whether it should be or not automatically

mask = bwareafilt(wshed(kMaskRefined), areaFilter);

%fit with circles instead of convex hulls
% [mask, ellipseFrame] = fitCirclesToNuclei(mask, kMask);

[mask, ellipseFrame] = fitCirclesToNuclei(mask, kMask);


%validate sizes. the ellipse masker handles
%very large objects poorly, to say the least.
if ~isempty(ellipseFrame)
    largeAxisIndex = ellipseFrame(:, 3) > max(size(image))...
        | ellipseFrame(3) > max(size(image));
    ellipseFrame(largeAxisIndex, :) = [];
end


end