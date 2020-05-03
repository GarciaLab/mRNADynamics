function [mask, ellipseFrame] = kSnakeCircles(image,...
    PixelSize_um, varargin)

%parameters i've found to be broadly applicable
sigmaK_um = .5;
mu = .1; %weight of length term for chen vese  algorithm. honestly don't know what this controls
min_rad_um = 1; % set min and max acceptable area for nucleus segmentation
max_rad_um = 6; %this needs to be 6um for nc12. 4um for nc14
nIterSnakes = 100;
maxAspectRatio = 4;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for k = 1:2:(numel(varargin)-1)
    if k ~= numel(varargin)
        eval([varargin{k} '=varargin{k+1};']);
    end
end


minArea_px = round(pi*(min_rad_um ./ PixelSize_um).^2);
maxArea_px = round(pi*(max_rad_um ./ PixelSize_um).^2);
sigmaK_px = sigmaK_um / PixelSize_um;
areaFilter = [minArea_px, maxArea_px];

%preprocessing denoise
image = wiener2(image);


%assuming the right class is 3 here. sometimes that might be wrong.
%it's hard to decide which one is right automatically.
%choose kMask does this. 
kLabel= imsegkmeans(single(imgaussfilt(image,sigmaK_px)),3);

kMask = kLabel == chooseKLabel(kLabel, areaFilter);

%sometimes snakes destroys blobs. if it does, it'd be nice to add back in
%regions from kMask.

kMask = imfill(kMask, 'holes');

kMask = bwareafilt(kMask, areaFilter); 

kMaskRefined= gather( chenvese( ...
    imgaussfilt( gpuArrayMaybe(image), sigmaK_px),...
    gpuArrayMaybe(kMask), nIterSnakes, mu, 'chan') );

%same issue here. sometimes the watershed is inverted. hard to pick
%whether it should be or not automatically

% kMaskRefined= imfill(kMaskRefined, 'holes');
mask = bwareafilt(wshed(kMaskRefined), areaFilter);

% mask = wshed(bwareafilt(kMaskRefined, areaFilter));

%morphologically clean up the mask
%removed because it introduces weird artifacts that connect
%far away regions
% mask = imfill(bwmorph(mask, 'bridge'), 'holes');

[mask, ellipseFrame] = fitCirclesToNuclei(mask, 'areaFilter', areaFilter, 'maxAspectRatio', maxAspectRatio);

mask = logical(mask); 


%validate sizes. the ellipse masker handles
%very large objects poorly

if ~isempty(ellipseFrame)
    largeAxisIndex = ellipseFrame(:, 3) > min(size(image))...
        | ellipseFrame(3) > min(size(image));
    ellipseFrame(largeAxisIndex, 3) = median(ellipseFrame(:, 3)); 
    ellipseFrame(largeAxisIndex, 4) = median(ellipseFrame(:, 4));
end

%validation
if ~isempty(ellipseFrame)
    assert( all(ellipseFrame(:, 5) <= 2*pi) );
end


figure(1); tiledlayout('flow');
nexttile; imagesc(image);
nexttile; imagesc(kMask);
nexttile; imagesc(kMaskRefined);
nexttile; imagesc(mask);
drawnow;


end