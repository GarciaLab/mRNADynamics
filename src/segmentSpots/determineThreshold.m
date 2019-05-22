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
    
    for i = 1:length(varargin)
        if strcmpi(varargin{i}, 'noSave')
            noSave = true;
            dogs = varargin{i+1};
        end
    end

    % loads information needed to loop through DOGs

    [SourcePath,ProcPath,DropboxFolder,MS2CodePath,PreProcPath]=...
        DetermineLocalFolders(Prefix);

    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo');
    zSize = FrameInfo(1).NumberSlices + 2;
    numFrames = length(FrameInfo);
    
    % says which z-slices and frames that we can scroll through
    available_zs = 2:3:(zSize - 2);
    available_frames = 1:4:numFrames;
    
    OutputFolder1=[ProcPath,filesep,Prefix,'_',filesep,'dogs',filesep];
    nameSuffix = ['_ch',iIndex(Channel,2)];
    
    % loops through DOGs to find brightest one
    % does this by choosing DOG with most pixels above brightest_iqr_test sds
    bestZ = 2;
    bestFrame = 1;
    bestVal = 0;
    max_val = 0;
    all_dogs = cell(numFrames, zSize - 2);
    dogDir = dir(OutputFolder1);
    dogDirNames = [dogDir.name];
    for frame = available_frames
        for z = available_zs
            dog = loadDog(z);
            all_dogs{frame, z - 1} = dog;
            non_zero_d = dog(dog ~= 0);
            val = iqr(non_zero_d(:))/2;
            median_val = median(non_zero_d(:));
            num_above = sum(non_zero_d(:) > (val * brightest_iqr_test + median_val));
            if num_above > bestVal
                bestVal = num_above;
                bestZ = z;
                bestFrame = frame;
            end
            max_val = max(max(dog(:)), max_val);
        end
    end
    % generates UI for picking threshold
    beiqrOG = all_dogs{bestFrame, bestZ - 1};
    non_zero_dog = beiqrOG(beiqrOG ~= 0);
    median_val = median(non_zero_dog(:));
    iqr_val = iqr(non_zero_dog(:))/2;
    dog_copy = all_dogs{bestFrame, bestZ - 1};
    min_val = min(non_zero_dog(:));
    thresh = min(median_val + default_iqr * iqr_val, max_val - 1);
    dog_copy(dog_copy < thresh) = 0;
    f = figure();
    uiAxes = axes(f);
    im = imshow(dog_copy, [], 'Parent', uiAxes);
    set(f, 'Position', [100, 100, 1000, 600])
    
    % Slider for changing the z slice
    zStep = 1.0 / (zSize - 2);
    zSlider = uicontrol('Style', 'slider', 'Min', 2, ...
        'Max', zSize - 1, 'Value', bestZ, 'SliderStep', [zStep, zStep], ...
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
    
    uiwait(f);
    
    function update_val(source, ~)
        
        zSlider.Value = round(zSlider.Value);
        frameSlider.Value = round(frameSlider.Value);
        bestZ = zSlider.Value;
        bestFrame = frameSlider.Value;
        if isempty(all_dogs{bestFrame, bestZ - 1})
            dog = loadDog(bestZ);
            all_dogs{bestFrame, bestZ - 1} = dog;
        end
        dog_copy = all_dogs{bestFrame, bestZ - 1};
        thresh = threshSlider.Value;
        iqr_above_median = (thresh - median_val) / iqr_val;
        dog_copy(dog_copy < thresh) = 0;
        im.CData = dog_copy;

        threshVal.String = ['threshold = ' num2str(thresh) ' which is ' ...
            num2str(iqr_above_median) ' sds above the median pixel value']; 
        zVal.String = ['Z-slice = ' num2str(bestZ)];
        frameVal.String = ['Frame = ' num2str(bestFrame)];
                 
        end
        

    function use_thresh(source, ~)
        uiresume(f);
        close(f);
    end

    function dog = loadDog(zInd)
        if ~noSave
            try
                dog_name = ['DOG_',Prefix,'_',iIndex(frame,3),'_z',iIndex(zInd,2),nameSuffix,'.tif'];
                dog = double(imread([OutputFolder1 dog_name]));
                
            catch
                dog_name = ['DOG_',Prefix,'_',iIndex(frame,3),'_z',iIndex(zInd,2),nameSuffix,'.mat'];
                load([OutputFolder1 dog_name], 'plane');
                dog = plane;
            end
            
        else
            dog = dogs(:, :, zInd, frame);
        end
    end
    
end

