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

cleanupObj = onCleanup(@myCleanupFun);

default_iqr = 6;
brightest_iqr_test = 8;
noSave = false;
dogs = [];
firstFrame = 1;
lastFrame = 0;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'noSave')
        noSave = true;
        dogs = varargin{i+1};
    elseif strcmpi(varargin{i}, 'numFrames') ||...
            strcmpi(varargin{i}, 'lastFrame')
        lastFrame = varargin{i+1};
    elseif strcmpi(varargin{i}, 'firstFrame')
        firstFrame = varargin{i+1};
    end
end

liveExperiment = LiveExperiment(Prefix);
ProcPath = liveExperiment.userProcFolder;
FrameInfo = getFrameInfo(liveExperiment);
zSize = FrameInfo(1).NumberSlices + 2;


if nargin < 2
    spotChannels = liveExperiment.spotChannels;
    if numel(spotChannels)==1
        Channel = spotChannels;
    else
        error('Please supply the spot channel as a second argument.');
    end
end

dogFolder=[ProcPath,filesep,Prefix,'_',filesep,'dogs',filesep];
dogDir = dir([dogFolder, filesep, '*_ch0', num2str(Channel), '.*']);
%get a cleaner directory with only preferred
%filetypes
isDogFolder = true;
dogDir = getImageFiletypes(dogDir, isDogFolder);

firstDogStack = imreadStack([dogFolder, dogDir(1).name]);

zPadded = zSize ~= size(firstDogStack, 3);


haveStacks = any(~contains(string({dogDir.name}), '_z'));

% determines which z-slices and
%frames that we can scroll through
if zPadded
    minZ = 2;
    maxZ = zSize - 1;
else
    minZ = 1;
    maxZ = zSize - 2;
end

skipZFactor = 1;
available_zs = 2:skipZFactor:zSize-1;


if lastFrame == 0
    lastFrame = numel(FrameInfo);
end

skipFrameFactor = 4;
available_frames = firstFrame:skipFrameFactor:lastFrame;

% loops through DOGs to find brightest one
% does this by choosing DOG
%with most pixels above brightest_iqr_test sds
bestZ = 2;
bestFrame = 1;
bestVal = 0;
max_val = 0;

disp('Loading dogs...')
all_dogs = cell(lastFrame, zSize - 2);
for frame = available_frames
    
    if haveStacks
        try dogStack = loadDogStack(frame);
        catch, continue; end
    end
    
    for z = available_zs
        
        if zPadded,  zInd = z;
        else, zInd = z-1; end
        
        if ~haveStacks
            dog = loadDogPlane(zInd, frame);
        else
            dog = dogStack(:, :, zInd);
        end
        
        all_dogs{frame, zInd} = dog;
        non_zero_d = dog(dog~=0);
        val = (1/2) * iqr( non_zero_d(:) );
        median_val = median( non_zero_d(:) );
        num_above = sum( non_zero_d(:) >...
            (val * brightest_iqr_test + median_val) );
        if num_above > bestVal
            bestVal = num_above;
            bestZ = zInd;
            bestFrame = frame;
        end
        max_val = max( max( dog(:) ), max_val);
        
    end
end
disp('Dogs loaded.');
%%
% get the image we'll display
bestDOG = all_dogs{bestFrame, bestZ};
non_zero_dog = bestDOG(bestDOG ~= 0);
median_val = median(non_zero_dog(:));
iqr_val = iqr(non_zero_dog(:))/2;
dog_copy = all_dogs{bestFrame, bestZ};
min_val = min(non_zero_dog(:));
thresh = min( median_val + (default_iqr * iqr_val), max_val);
dog_copy = dog_copy > thresh;

nSpotsEstimate = length(regionprops(dog_copy));

%%
% generates UI for picking threshold
fig = figure();
uiAx = axes(fig);
im = imshow(dog_copy, [], 'Parent', uiAx);

% lays out user interface
screen_size = get(0, 'screensize');
dim = [screen_size(3) * 0.6, screen_size(4) * 0.75];
dimVec = [dim(1), dim(2), dim(1), dim(2)]; %to easily normalize units

set(fig, 'Position', [100, 100, 1000, 600])

% Slider for changing the z slice
zStep = 1.0 / (zSize - 2);

zSlider = uicontrol('Style', 'slider', 'Min', minZ, ...
    'Max', maxZ, 'Value', bestZ, 'SliderStep', [zStep, zStep], ...
    'Position', [250 50 500 20], 'Callback', @update_val);

zVal = uicontrol('Style','text',...
    'Position',[750 50 100 20],...
    'String',['Z-slice = ' num2str(bestZ)]);

% Slider for changing the frame
fStep = 1.0 / lastFrame;
frameSlider = uicontrol('Style', 'slider', 'Min', 1, ...
    'Max', lastFrame, 'Value', bestFrame, 'SliderStep', [fStep, fStep], ...
    'Position', [60 50 20 500], 'Callback', @update_val);

frameVal = uicontrol('Style','text',...
    'Position',[60 555 100 20],...
    'String',['Frame = ' num2str(bestFrame)]);

% Slider for adjusting the threshold
threshSlider = uicontrol('Style', 'slider',...
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


%% display contrast stuff
maxPos =  [.75, .3, .2, .03];
maxLabelPos =[maxPos(1), maxPos(2)+.025,...
    maxPos(3)/2, maxPos(4)+.02];
%make sure the display range sliders and labels
%move and resize as a group
minPos =  [maxPos(1), maxPos(2) + .1, maxPos(3), maxPos(4)];
minLabelPos = [minPos(1), minPos(2)+.025,...
    maxLabelPos(3), maxLabelPos(4)];

display_range = double([min_val, max_val]);

max_label = uicontrol(fig,  'Style','text', 'String', 'max display value',...
    'Units', 'normalized', 'Position', maxLabelPos);
max_slider = uicontrol(fig, 'Style', 'slider', 'Min', display_range(1),...
    'Max', display_range(2), 'Value', display_range(2), ...
    'Units', 'normalized', 'Position', maxPos,...
    'Callback', @update_val);

min_label = uicontrol(fig,  'Style','text', 'String', 'min display value', ...
    'Units', 'normalized', 'Position', ...
    minLabelPos);
min_slider = uicontrol(fig, 'Style', 'slider', 'Min', display_range(1),...
    'Max', display_range(2), 'Value', display_range(1), ...
    'Units', 'normalized', 'Position', minPos,...
    'Callback', @update_val);
%%

nSpotsLabelPos =  [minPos(1),...
    minLabelPos(2)+.1, minLabelPos(3), minLabelPos(4)];

nSpotsLabel = uicontrol(fig, 'Style', 'text','String',...
    ['Estimated number of spots: ', num2str(nSpotsEstimate) ],...
    'Units', 'normalized', 'Position', nSpotsLabelPos);

% 
% nSpotsAxes = axes(fig, 'Units', 'normalized', 'Position', [nSpotsLabelPos(1),...
%     nSpotsLabelPos(2)+.075,  .2, .2]);

cboxPos = [825 100 200 20];

%checkbox for displaying unthresholded image
chk = uicontrol('Style', 'checkbox', 'String', 'Display unthresholded image',...
    'Position',cboxPos,...
    'Callback', @update_val);
%%
uiwait(fig);

    function update_val(source, ~)
        
        zSlider.Value = round(zSlider.Value);
        frameSlider.Value = round(frameSlider.Value);
        if frameSlider.Value <= 0
            frameSlider.Value = 1;
        end
        bestZ = zSlider.Value;
        bestFrame = frameSlider.Value;
        if isempty(all_dogs{bestFrame, bestZ})
            try
            dog = loadDogPlane(bestZ,bestFrame);
            catch, return; 
            end
            
            all_dogs{bestFrame, bestZ} = dog;
        end
        dog_copy = all_dogs{bestFrame, bestZ};
        thresh = threshSlider.Value;
        iqr_above_median = (thresh - median_val) / iqr_val;
        
        maxDisplayIntensity = round(max_slider.Value);
        minDisplayIntensity = round(min_slider.Value);
        
        if minDisplayIntensity >= maxDisplayIntensity
            maxDisplayIntensity = max(max(dog_copy));
            minDisplayIntensity = median(median(dog_copy));
        end
        
        dogbw = dog_copy > thresh;
        
        %         %%
        %         nSpotsF = zeros(1, length(available_frames));
        %         for f = available_frames
        %             nSpotsZ = zeros(1, length(available_zs));
        %             for iz = available_zs
        %                 if zPadded
        %                     zIndex = iz;
        %                 else
        %                     zIndex = iz-1;
        %                 end
        %                 if isempty(all_dogs{f, zIndex})
        %                     dog = loadDogPlane(f, zIndex);
        %                     all_dogs{f, zIndex} = dog;
        %                 end
        %                 nSpotsZ(zIndex) = length(regionprops(all_dogs{f, zIndex} > thresh));
        %
        %             end
        %             nSpotsF(f)= max(nSpotsZ);
        %         end
        %         nSpotsEstimate = length(regionprops(dogbw));
        %         nSpotsLabel.String = ['Estimated number of spots: ', num2str(nSpotsEstimate) ];
        %
        %         if ~isempty(nSpotsAxes.Children)
        %             nSpotsAxes.Children(1).YData = nSpotsF(available_frames);
        %             nSpotsAxes.Children(1).XData = available_frames;
        %         else
        %             line(nSpotsAxes, available_frames, nSpotsF(available_frames))
        %         end
        %
        %%
        drawnow;
        
        if ~chk.Value
            dog_copy = dogbw;
            im.CData = dog_copy;
            uiAx.CLim = [0, 1];
            drawnow;
        else
            im.CData = dog_copy;
            im.CDataMapping = 'scaled';
            uiAx.CLim = [minDisplayIntensity, maxDisplayIntensity];
            drawnow;
        end
        
        threshVal.String = ['threshold = ' num2str(thresh) ' which is ' ...
            num2str(iqr_above_median) ' IQRs above the median pixel value'];
        zVal.String = ['Z-slice = ' num2str(bestZ)];
        frameVal.String = ['Frame = ' num2str(bestFrame)];
        
    end


    function use_thresh(source, ~)
        uiresume(fig);
        close(fig);
    end

    function dog = loadDogPlane(zInd, frame)
        
        if haveStacks
            
            filenameIndex = contains(string({dogDir.name}), iIndex(frame, 3));
            file = [dogFolder, dogDir(filenameIndex).name];
            dogStack = imreadStack(file);
            dog = dogStack(:, :, zInd);
            
        else
            
            filenameIndex = contains(string({dogDir.name}), iIndex(frame, 3)) &...
            contains(string(dogDir.name), iIndex(zInd, 2));
            file = [dogFolder, dogDir(filenameIndex).name];
            
            if contains(file, '.tif')
                dog = double(imread(file));
            elseif contains(file, '.mat')
                load(file, 'plane');
                try dog = plane; 
                catch, load(file, 'dog'); end
            end
            
        end
        
    end

    function dogStack = loadDogStack(frame)
        
            filenameIndex = contains(string({dogDir.name}), iIndex(frame, 3));
            file = [dogFolder, dogDir(filenameIndex).name];
            dogStack = imreadStack(file);
            
    end 

end
