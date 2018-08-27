% determineThreshhold(Prefix)
%
% DESCRIPTION
% Tool for finding threshhold for segmentSpots. Generates a UI and starts
% by picking a threshhold 4 standard deviations about the mean in the
% brightest image. You can change this threshhold with a slider. Everything
% below the threshhold is zeroes out so you can see what would be picked up
% by segmentSpots
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% Channel: Channel to get threshhold for
%
% Author (contact): Sean Medin (smedin@berkeley.edu)
% Created: 8/23/2018
%
% Documented by: Sean Medin (smedin@berkeley.edu)


function [thresh] = determineThreshhold(Prefix, Channel)

    default_std = 5;
    num_frames_to_check = 50;

    % loads information needed to loop through DOGs
    [~,~,~,~,~,~,~,ExperimentType, Channel1, Channel2,~] =...
    readMovieDatabase(Prefix);

    [SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
        DetermineLocalFolders(Prefix);

    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo');
    zSize = FrameInfo(1).NumberSlices + 2;
    numFrames = length(FrameInfo);
    OutputFolder1=[FISHPath,filesep,Prefix,'_',filesep,'dogs',filesep];
    nameSuffix = ['_ch',iIndex(Channel,2)];
    
    % loops through DOGs to find brightest one
    bestDOG = 0;
    bestVal = 0;
    num_skip = ceil(numFrames / num_frames_to_check); % looks at 50 frames only
    for frame = 1:num_skip:numFrames
        for z = 2:(zSize - 1)
            
            dog_name = ['DOG_',Prefix,'_',iIndex(frame,3),'_z',iIndex(z,2),nameSuffix,'.tif'];
            dog = double(imread([OutputFolder1 dog_name]));
            val = sum(sum(dog));
            if val > bestVal
                bestVal = val;
                bestDOG = dog;
            end
        end
    end
    % generates UI for picking threshhold
    mean_val = mean(bestDOG(:));
    std_val = std(bestDOG(:));
    dog_copy = bestDOG;
    thresh = mean_val + default_std * std_val;
    dog_copy(dog_copy < thresh) = 0;
    f = figure();
    imshow(dog_copy, []);
    set(gcf, 'Position', [100, 100, 1000, 600])
    threshSlider = uicontrol('Style', 'slider',...
        'Min',1,'Max',max(bestDOG(:)),'Value',thresh,...
        'Position', [200 20 500 20],...
        'Callback', @update_val); 
    
    threshVal = uicontrol('Style','text',...
        'Position',[200 45 500 20],...
        'String',['Threshhold = ' num2str(thresh) ' which is ' ...
        num2str(default_std) ' sds above the mean pixel value']);
    
    % Create push button
    btn = uicontrol('Style', 'pushbutton', 'String', 'Use Threshhold',...
        'Position', [700 20 250 20],...
        'Callback', @use_thresh); 
    
    uiwait(gcf);
    
    function update_val(source, ~)
        dog_copy = bestDOG;
        thresh = source.Value;
        std_above_mean = (thresh - mean_val) / std_val;
        dog_copy(dog_copy < thresh) = 0;
        imshow(dog_copy, []);
        set(gcf, 'Position', [100, 100, 1000, 600])
        threshSlider = uicontrol('Style', 'slider',...
            'Min',1,'Max',max(bestDOG(:)),'Value',thresh,...
            'Position', [200 20 500 20],...
            'Callback', @update_val); 

        threshVal = uicontrol('Style','text',...
            'Position',[200 45 500 20],...
            'String',['Threshhold = ' num2str(thresh) ' which is ' ...
            num2str(std_above_mean) ' sds above the mean pixel value']); 
        
        btn = uicontrol('Style', 'pushbutton', 'String', 'Use Threshhold',...
        'Position', [700 20 250 20],...
        'Callback', @use_thresh);
        
    end

    function use_thresh(source, ~)
        uiresume(gcf);
        close(f);
    end
    
end

