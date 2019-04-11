function [ xy, varargout ] = findNuclei(FrameInfo, names, frameNumber, nucleusDiameter, embryoMask, varargin )
%FINDNUCLEI This function finds the position of nuclei within a frame which
%is used as an initial input for the tracking algorithm.
%   Detection is made by detecting the local maxima in a
%   Laplacian-of-Gaussian filtered image.
%

%% Initializing variables
% Load parameters
numberOfFrames = numel(names);
LoGratio = getDefaultParameters(FrameInfo,'LoGratio');
space_resolution = getDefaultParameters(FrameInfo,'space resolution');
localMaximumRadius = LoGratio*nucleusDiameter/space_resolution;
LoGradius = nucleusDiameter/space_resolution*LoGratio;
edgeClearance = getDefaultParameters(FrameInfo,'edge clearance')*nucleusDiameter/space_resolution;

if ~exist('embryoMask','var') || isempty(embryoMask)
    embryoMask = true(size(imread(names{frameNumber})));
end

if nargin > 5
    targetNumber = varargin{1}; % coarse estimate of the number of nuclei that should be found.
end

% Create the mask used to detect local maxima.
localMaxMask = fspecial('disk',localMaximumRadius);
localMaxMask = im2bw(mat2gray(localMaxMask),graythresh(mat2gray(localMaxMask)));
localMaxMask(round(length(localMaxMask)/2),round(length(localMaxMask)/2))  = 0;


%% Main body

% Load image
img = double(imread(names{frameNumber}));

% Filter the image
%filteredImg = imfilter(img,-fspecial('log',round(10*LoGradius),LoGradius),'symmetric');
% filteredImg = imresize(fourierFilterWithSymmetricBoundaryConditions(imresize(img,.5),-fspecial('log',round(10*LoGradius),LoGradius)), 2);
filteredImg = fourierFilterWithSymmetricBoundaryConditions(img,-fspecial('log',round(10*LoGradius),LoGradius));

% Find local maxima
% maxima = (filteredImg > imresize(imdilate(imresize(filteredImg, .5),imresize(localMaxMask, .5)),2) ) & embryoMask;
maxima = (filteredImg > imdilate(filteredImg,localMaxMask) ) & embryoMask;

% Smooth the image to get more robust maxima values:
G = imfilter(filteredImg,fspecial('disk',3),'symmetric');
ind = find(maxima>0);
[xm,ym] = ind2sub(size(maxima),ind);

%Check whether anything was found in this image. If nothing is found,
%this is usually the result of an empty frame
if isempty(xm)|isempty(ym)
   error(['No nuclei found in frame ',num2str(frameNumber),'. Check that that frame is not blank.']) 
end

try
    [v,dummy] = voronoin([xm ym]);
catch
    [v,dummy] = voronoin([xm ym],{'Qbb','Qz'});
end

ind_v = v(:,1) > 0.5 & v(:,2) > 0.5 & v(:,1) < size(maxima,1) & v(:,2) < size(maxima,2);
background_ind = sub2ind(size(maxima),round(v(ind_v,1)),round(v(ind_v,2)));
sample_ind = [ind(:) ; background_ind(:)];
nuc1 = G(sample_ind);

thresh = graythresh(mat2gray(nuc1));
indNuclei = im2bw(mat2gray(nuc1),thresh);

%% Check the segmentation if a target number was provided.
if exist('targetNumber','var') && numel(targetNumber) == 1 && isnumeric(targetNumber) && abs(sum(indNuclei(1:numel(nuc1)))-targetNumber)/targetNumber > 0.25
    % Enforce a certain number of nuclei on the current frame, plus or
    % minus 'perc' percent
    perc = 25;
    thresh = graythresh(mat2gray(nuc1));
    indNuclei = im2bw(mat2gray(nuc1),thresh);
    % Test whether the number of nuclei found is within the given range. If
    % not, then chose the threshold to get as close as possible to the 
    % target number.
    if abs(sum(indNuclei(1:numel(nuc1)))-targetNumber)/targetNumber > perc/100;
        indNucleiTmp = nan(21,numel([nuc1;nuc_prev;nuc_next]));
        T = 0:.005:1;
        % test a range of threshold and pick the one that returns the
        % number of nuclei the closest to the target number.
        for j = 1:numel(T)
            indNucleiTmp(j,:) = im2bw(mat2gray([nuc1;nuc_prev;nuc_next]),T(j));
        end
        indNucleiTmp(:,(numel(nuc1)+1):end) = [];
        [dummy, IND] =  min(abs(sum(indNucleiTmp,2)-targetNumber));
        indNuclei = logical(indNucleiTmp(IND,:));
    end
end
indNuclei = indNuclei(1:numel(ind));

% Get rid of all nuclei below thresh
ind = ind(indNuclei);

% Convert into coordinates
[x,y] = ind2sub(size(img),ind);

% Get rid of all nuclei too close to edges
%indNuclei = x >= 1+edgeClearance & x <= size(img,1)-edgeClearance & y >= 1+edgeClearance & y <= size(img,2)-edgeClearance;
%
%x = x(indNuclei);
%y = y(indNuclei);


xy = [x y];

if nargout > 1
    varargout{1} = G(ind); % Those values are the intensities of the maxima that could be used to sort them.
end
end

