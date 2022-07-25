function normalizeForBleaching(Prefix, channel)
%
% DESCRIPTION
% This function takes all the 3D (xyz) TIFF frames from a single movie and 
% normalizes their pixel values to the mean-above-background pixel
% intensity
% Note: This function only works with the new 3D TIFF stacks that are
%       output by the current version of the pipeline
%
% MT: I wrote this function to create movies that are easier for Weka to
%     classify using a single classifier for the whole movie, even if
%     bleaching is severe.
%
% INPUT ARGUMENTS
% prefix:
%
% channel: 
% 
% OPTIONS
% N/A
%
% OUTPUT
% normalizedFolder: path to the folder where the normalized movie frames
%                   (saved as 3D, xyz TIFF stacks, same as the input movie)
%                   are saved
% TIFF files containing a bleaching corrected movie
%
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 08/12/2020
% Last Updated: N/A
%

% channel = 1;
% Get needed info for this experiment
liveExperiment = LiveExperiment(Prefix);
preProcFolder = liveExperiment.preFolder;

ch = num2str(channel);
rawImDir{channel} = dir([preProcFolder '*ch0', ch, '.tif']);
normImWritePath = [preProcFolder 'normalizedImages\'];
mkdir(normImWritePath);

%% Get the mean pixel value above background from the first frame of the movie
firstImPath = [rawImDir{1,channel}(1).folder, filesep, rawImDir{1,channel}(1).name];
firstImStack = loadTiffStack(firstImPath);
zDim = size(firstImStack,3);
% Get the mean intensity above background for each z slice
for z = 1:zDim
    currZSlice = firstImStack(:,:,z);
    meanIntensity(1,z) = mean(currZSlice(:));
    meanAboveBackground(1,z) = mean(currZSlice(currZSlice > mean(currZSlice(:)))); % get background
end
% for simplicity (and because the z-stack probably was shifted during data
% taking), take the average of all the z-slice means
meanAboveBackground0 = nanmean(meanAboveBackground(1,:));

% Copy first frame stack, unmodified, to the normalized movie folder
copyfile(firstImPath, normImWritePath);

% Save max projected image for easy comparison of consecutive frames
firstImZMaxProj = max(firstImStack,[],3);
% firstImZSumProj = sum(firstImStack,3);
% imshow(firstImZSumProj,[])
nameSuffix = ['_ch',iIndex(channel,2)];
imageMaxName = [Prefix, '_', iIndex(1,3), '_normMax', nameSuffix, '.tif'];
mkdir([normImWritePath, 'normMaxProj', filesep]);
imwrite(uint16(firstImZMaxProj), [normImWritePath, 'normMaxProj', filesep, imageMaxName]);


%% Normalize the remaining frames in the movie to the first frame
nFrames = numel(rawImDir{channel});
for currFrame = 2:nFrames
    %Get the raw image
    currImPath = [rawImDir{1,channel}(currFrame).folder, filesep, rawImDir{1,channel}(currFrame).name];
    currImStack = loadTiffStack(currImPath); %this is the xDim x yDim x zDim image matrix
    
    nameSuffix = ['_ch',iIndex(channel,2)];
    imageName = [Prefix, '_', iIndex(currFrame,3), '_norm', nameSuffix, '.tif'];
    
    % Indivually normalize each z slice
    for z = 1:zDim
        % 1st slice is padding/blank, no need to normalize
        if z == 1
            firstZSlice = currImStack(:,:,1);
            normImStack(:,:,1) = firstZSlice;    
            imwrite(uint16(firstZSlice), [normImWritePath, filesep, imageName]);
        % last slice is padding/blank, no need to normalize
        elseif z == zDim
            lastZSlice = currImStack(:,:,zDim);
            normImStack(:,:,zDim) = lastZSlice; 
            imwrite(uint16(lastZSlice),[normImWritePath, filesep, imageName], 'WriteMode', 'append');
        else
            currZSlice = currImStack(:,:,z);
            % find the mean intensity above background
            meanIntensity(currFrame,z) = mean(currZSlice(:));
            meanAboveBackground(currFrame,z) = mean(currZSlice(currZSlice > meanIntensity(currFrame,z)));
            % normalize the intensity of this slice by the 1st frame's mean
            normRatio(currFrame,z) = meanAboveBackground0 / meanAboveBackground(currFrame,z);
            normZSlice = currZSlice .* normRatio(currFrame,z);
            normImStack(:,:,z) = normZSlice;
            % append normalized frame to tiff stack
            imwrite(uint16(normZSlice),[normImWritePath, filesep, imageName], 'WriteMode', 'append');
        end
    end
  
    imZSumProj = sum(normImStack,3);
    imZMaxProj = max(normImStack,[],3);
    zSumProjArray(:,:,currFrame) = imZSumProj;
    zMaxProjArray(:,:,currFrame) = imZMaxProj;
    
    % Display and save max projected, normalized frame
%     imshow(imZMaxProj,[])
%     drawnow
    imageMaxName = [Prefix, '_', iIndex(currFrame,3), '_normMax', nameSuffix, '.tif'];
    imwrite(uint16(imZMaxProj),[normImWritePath, 'normMaxProj', filesep, imageMaxName], 'WriteMode', 'append');
end