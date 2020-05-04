function [ xy, varargout ] = findNuclei(FrameInfo, hisMat,...
    frameNumber, nucleusDiameter, embryoMask, varargin )
%FINDNUCLEI This function finds the position of nuclei within a frame which
%is used as an initial input for the tracking algorithm.
%   Detection is made by detecting the local maxima in a
%   Laplacian-of-Gaussian filtered image.
%

%% Initializing variables
% Load parameters

LoGratio = getDefaultParameters(FrameInfo,'LoGratio');
space_resolution = getDefaultParameters(FrameInfo,'space resolution');
localMaximumRadius = LoGratio*nucleusDiameter/space_resolution;
LoGradius = nucleusDiameter/space_resolution*LoGratio;
edgeClearance = getDefaultParameters(FrameInfo,'edge clearance')*nucleusDiameter/space_resolution;
%% 

% Added by NL and GM on 11/23/2019
% Edited by GM on 1/7/2019
xDim = FrameInfo(1).PixelsPerLine * FrameInfo(1).PixelSize;
yDim = FrameInfo(1).LinesPerFrame * FrameInfo(1).PixelSize;
if yDim > 150 && xDim > 150
    I = double(hisMat(:, :, frameNumber));
    pixelvalues = unique(I(:));
    thresh = pixelvalues(2);
    f_sigma = round(nucleusDiameter / FrameInfo(1).PixelSize);
    bwfill = imfill(I>thresh,'holes');
    I_inside = bwselect(bwfill,round(size(bwfill,2)/2), round(size(bwfill, 1)/2));
    I_inside=uint16((2^16-1)*I_inside);
    I_blurred = imfilter(I_inside,...
         fspecial('gaussian',2*f_sigma,f_sigma),'symmetric','conv');
    level = graythresh(I_blurred);
    embryoMask = im2bw(I_blurred,level);
else    
    if ~exist('embryoMask','var') || isempty(embryoMask)
        embryoMask = true(size(hisMat(:, :, frameNumber), 1));
    end
end

%% 

if nargin > 5
    targetNumber = varargin{1}; % coarse estimate of the number of nuclei that should be found.
end
    

% Create the mask used to detect local maxima.
localMaxMask = fspecial('disk',localMaximumRadius);
localMaxMask = imbinarize(mat2gray(localMaxMask),graythresh(mat2gray(localMaxMask)));
localMaxMask(round(length(localMaxMask)/2),round(length(localMaxMask)/2))  = 0;


%% Main body

% Load image
img = double(hisMat(:, :, frameNumber));

% Filter the image
filteredImg = fourierFilterWithSymmetricBoundaryConditions(img,-fspecial('log',round(10*LoGradius),LoGradius));

% Find local maxima
maxima = (filteredImg > imdilate(filteredImg,localMaxMask) ) & embryoMask;

% Smooth the image to get more robust maxima values:
G = imfilter(filteredImg,fspecial('disk',3),'symmetric');
ind = find(maxima>0);
[xm,ym] = ind2sub(size(maxima),ind);

%Check whether anything was found in this image. If nothing is found,
%this is usually the result of an empty frame
if isempty(xm) || isempty(ym)
   error(['No nuclei found in frame ',num2str(frameNumber),'. Check that that frame is not blank.']) 
end

try
    [v,~] = voronoin([xm ym]);
catch
    [v,~] = voronoin([xm ym],{'Qbb','Qz'});
end

ind_v = v(:,1) > 0.5 & v(:,2) > 0.5 & v(:,1) < size(maxima,1) & v(:,2) < size(maxima,2);
background_ind = sub2ind(size(maxima),round(v(ind_v,1)),round(v(ind_v,2)));
sample_ind = [ind(:) ; background_ind(:)];
nuc1 = G(sample_ind);

thresh = graythresh(mat2gray(nuc1));
indNuclei = imbinarize(mat2gray(nuc1),thresh);

%% Check the segmentation if a target number was provided.
if exist('targetNumber', 'var')
    indNuclei = useTargetNumber(targetNumber, nuc1, indNuclei);
end
%%

indNuclei = indNuclei(1:numel(ind));

% Get rid of all nuclei below thresh
ind = ind(indNuclei);

% Convert into coordinates
[x,y] = ind2sub(size(img),ind);

xy = [x y];

if nargout > 1
    varargout{1} = G(ind); % Those values are the intensities of the maxima that could be used to sort them.
end

end

