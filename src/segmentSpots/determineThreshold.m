% determineThreshold(Prefix,Channel)
%
% DESCRIPTION
% Tool for finding threshold for segmentSpots. Generates a UI and starts
% by picking a threshold 4 standard deviations about the median in the
% brightest image. You can change this threshold with a slider. Everything
% below the threshold is zeroes out so you can see what would be picked up
% by segmentSpots
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% Channel: Channel to get threshold for
%
% Author (contact): Sean Medin (smedin@berkeley.edu)
% Created: 8/23/2018
%
% Documented by: Sean Medin (smedin@berkeley.edu)


function [thresh] = determineThreshold(Prefix, Channel, varargin)

default_iqr = 6;
brightest_iqr_test = 8;
noSave = false;
dogs = [];

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'noSave')
        noSave = true;
        dogs = varargin{i+1};
    elseif strcmpi(varargin{i}, 'numFrames')
        numFrames = varargin{i+1};
    end
end

% loads information needed to loop through DoGs

[~,ProcPath,DropboxFolder,~,~]=...
    DetermineLocalFolders(Prefix);

load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo');
zSize = FrameInfo(1).NumberSlices + 2;

OutputFolder1=[ProcPath,filesep,Prefix,'_',filesep,'dogs',filesep];

DogOutputFolder= OutputFolder1;

if nargin < 2
    spotChannels = thisExperiment.getSpotChannels; 
    if numel(spotChannels)==1 
        Channel = spotChannels;
    else
        error('Please supply the spot channel as a second argument.');
    end
end
dogDir = dir([OutputFolder1, filesep, '*_ch0', num2str(Channel), '.*']);

loadAsStacks = ~contains(dogDir(1).name, '_z');
Weka = startsWith(dogDir(1).name, 'prob');
loadAsMat = endsWith(dogDir(1).name, '.mat');
% dogStr = 'dogStack_';

if Weka
    dogStr = 'prob';
elseif loadAsStacks
    dogStr = 'dogStack_';
else
    dogStr = 'DOG_';
end


nameSuffix = ['_ch',iIndex(Channel,2)];


firstDogStackFile = [OutputFolder1, filesep, dogStr, Prefix, '_', iIndex(1, 3),...
    nameSuffix];
if loadAsMat
    load([firstDogStackFile, '.mat'], 'dogStack');
else
    dogStack = imreadStack([firstDogStackFile, '.tif']);
end

zPadded = zSize ~= size(dogStack, 3);

% says which z-slices and frames that we can scroll through
if zPadded
    minZ = 2;
    maxZ = zSize - 1;
else
    minZ = 1;
    maxZ = zSize - 2;
end
available_zs = 2:3:zSize-1;


numFrames = numel(dogDir);

if numFrames == 0
    numFrames = numel(FrameInfo);
end

available_frames = 1:4:numFrames;

% loops through DOGs to find brightest one
% does this by choosing DOG with most pixels above brightest_iqr_test sds
bestZ = 2;
bestFrame = 1;
bestVal = 0;
max_val = 0;
all_dogs = cell(numFrames, zSize - 2);
for frame = available_frames
    for z = available_zs
        if zPadded
            zInd = z;
        else
            zInd = z-1;
        end
        dog = loadDog(zInd, frame);
        all_dogs{frame, zInd} = dog;
        non_zero_d = dog(dog ~= 0);
        val = iqr(non_zero_d(:))/2;
        median_val = median(non_zero_d(:));
        num_above = sum(non_zero_d(:) > (val * brightest_iqr_test + median_val));
        if num_above > bestVal
            bestVal = num_above;
            bestZ = zInd;
            bestFrame = frame;
        end
        max_val = max(max(dog(:)), max_val);
    end
end
% generates UI for picking threshold
bestDOG = all_dogs{bestFrame, bestZ};
non_zero_dog = bestDOG(bestDOG ~= 0);
median_val = median(non_zero_dog(:));
iqr_val = iqr(non_zero_dog(:))/2;
dog_copy = all_dogs{bestFrame, bestZ};
min_val = min(non_zero_dog(:));
thresh = min(median_val + default_iqr * iqr_val, max_val - 1);
dog_copy(dog_copy < thresh) = 0;


f = figure();
uiAxes = axes(f);
im = imshow(dog_copy, [], 'Parent', uiAxes);
set(f, 'Position', [100, 100, 1000, 600])

% Slider for changing the z slice
zStep = 1.0 / (zSize - 2);

zSlider = uicontrol('Style', 'slider', 'Min', minZ, ...
    'Max', maxZ, 'Value', bestZ, 'SliderStep', [zStep, zStep], ...
    'Position', [250 50 500 20], 'Callback', @update_val);

zVal = uicontrol('Style','text',...
    'Position',[750 50 100 20],...
    'String',['Z-slice = ' num2str(bestZ)]);

% Slider for changing the frame
fStep = 1.0 / numFrames;
frameSlider = uicontrol('Style', 'slider', 'Min', 1, ...
    'Max', numFrames, 'Value', bestFrame, 'SliderStep', [fStep, fStep], ...
    'Position', [60 50 20 500], 'Callback', @update_val);

frameVal = uicontrol('Style','text',...
    'Position',[60 555 100 20],...
    'String',['Frame = ' num2str(bestFrame)]);

% Slider for adjusting the threshold
threshSlider = uicontrol('Style', 'slider',... 332q
    'Min',min_val,'Max',max_val,'Value',thresh,...
    'Position', [250 5 500 20],...
    'Callback', @update_val);

threshVal = uicontrol('Style','text',...
    'Position',[250 25 500 20],...
    'String',['threshold = ' num2str(thresh) ' which is ' ...
    num2str(default_iqr) ' sds above the median pixel value']);

% Button to press once threshold is found
btn = uicontrol('Style', 'pushbutton', 'String', 'Use threshold',...
    'Position', [750 5 250 20],...
    'Callback', @use_thresh);

%checkbox for displaying unthresholded image
chk = uicontrol('Style', 'checkbox', 'String', 'Display unthresholded image',...
    'Position', [825 100 200 20],...
    'Callback', @update_val);

uiwait(f);

    function update_val(source, ~)
        
        zSlider.Value = round(zSlider.Value);
        frameSlider.Value = round(frameSlider.Value);
        if frameSlider.Value <= 0
            frameSlider.Value = 1;
        end
        bestZ = zSlider.Value;
        bestFrame = frameSlider.Value;
        if isempty(all_dogs{bestFrame, bestZ})
            dog = loadDog(bestZ,bestFrame);
            all_dogs{bestFrame, bestZ} = dog;
        end
        dog_copy = all_dogs{bestFrame, bestZ};
        thresh = threshSlider.Value;
        iqr_above_median = (thresh - median_val) / iqr_val;
        if ~chk.Value
            dog_copy(dog_copy < thresh) = 0;
            im.CData = dog_copy;
        else
            im.CData = dog_copy;
            im.CDataMapping = 'scaled';
            uiAxes.CLim = [median(median(dog_copy)), max(max(dog_copy))];
            %             im = imagescUpdate(uiAxes,dog_copy,[median(median(dog_copy)), max(max(dog_copy))]);
        end
        
        threshVal.String = ['threshold = ' num2str(thresh) ' which is ' ...
            num2str(iqr_above_median) ' sds above the median pixel value'];
        zVal.String = ['Z-slice = ' num2str(bestZ)];
        frameVal.String = ['Frame = ' num2str(bestFrame)];
        
    end


    function use_thresh(source, ~)
        uiresume(f);
        close(f);
    end

    function dog = loadDog(zInd, frame)
        
        dogStackFile = [OutputFolder1, filesep, dogStr, Prefix, '_', iIndex(frame, 3),...
            nameSuffix];
        if loadAsMat
            load([dogStackFile, '.mat'], 'dogStack');
        else
            dogStack = imreadStack([dogStackFile, '.tif']);
%         dog_name = [dogProb,Prefix,'_',iIndex(frame,3),'_z',iIndex(zInd,2),nameSuffix,saveType];
        
%         if strcmpi(saveType, '.tif')
%             dog = double(imread([OutputFolder1 dog_name]));
%         elseif strcmpi(saveType, '.mat')
%             load([OutputFolder1 dog_name], 'plane', 'dog');
%             try dog = plane; end
%         elseif strcmpi(saveType, 'none')
%             dog = dogs(:, :, zInd, frame);
        end
        
        dog = dogStack(:, :, zInd);
        
        %         dog = double(squeeze(dogMat(:, :, zInd, frame)));

    end

end

