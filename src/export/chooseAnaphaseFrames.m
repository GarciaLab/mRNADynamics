% DESCRIPTION
% This function generates a GUI to explore different ways of generating a
% nuclear channel. You can invert a channel and you can combine multiple
% together. You can also change the projectiontype that is used. The final
% combination of channels and projection type used when pressing 'Save
% Channel Selection' will be what will be used when forming the
% nuclear/histone channel. Note that this does not change your
% MovieDatabase entry--just the nuclear/histone images generated. At the
% moment this function should only be called downstream of
% exportDataForLivemRNA and is only usable with Leica data (and the
% 'nuclearGUI' option must be entered in exportDataForLivemRNA)

function [anaphaseFrames, projectionChannels, ProjectionType, hisMat] =...
    chooseAnaphaseFrames(Prefix, varargin)

cleanupObj = onCleanup(@myCleanupFun);
warning('off', 'MATLAB:ui:Slider:fixedHeight')
warning('off', 'MATLAB:audiovideo:audioplayer:noAudioOutputDevice');
clear getHisMat;


skip_factor = 1; % Only uses 1/skip_factor frames

% default custom projection parameters
max_custom = 1; % highest histone channel slice used
min_custom = 5; % lowest histone channel slice used


ProjectionType = 'midsumprojection';
load('ReferenceHist.mat', 'ReferenceHist');

returnHisMat = nargout > 2;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for arg = 1:2:(numel(varargin)-1)
    if arg ~= numel(varargin)
        eval([varargin{arg} '=varargin{arg+1};']);
    end
end


liveExperiment = LiveExperiment(Prefix);

if ~isempty(Prefix) 
   
      
    Channel1 = liveExperiment.Channel1;
    Channel2 = liveExperiment.Channel2;
    Channel3 = liveExperiment.Channel3;
    
    projectionChannels = {Channel1, Channel2, Channel3};
    
    DropboxFolder = liveExperiment.userResultsFolder;
    PreProcPath = liveExperiment.userPreFolder;
    
    movieMat = getMovieMat(liveExperiment);
    
    projectionTypeFile = [DropboxFolder,filesep,Prefix,filesep, 'ProjectionType.mat'];
    
    if exist(projectionTypeFile, 'file')
        load(projectionTypeFile, 'ProjectionType')
    else
        ProjectionType = 'midsumprojection';
    end
    
    
    
    [anaphaseFrames, anaphaseFile] = retrieveAnaphaseFrames(Prefix);
    anaphaseFramesInitial = anaphaseFrames;
    
    channelsFile = [DropboxFolder,filesep,Prefix,filesep, 'Channels.mat'];

    
    isUnhealthyFile = [DropboxFolder,filesep,Prefix,filesep, 'isUnhealthy.mat'];
    if exist(isUnhealthyFile, 'file')
        load(isUnhealthyFile, 'isUnhealthy');
    end
    
    audioFile = 'X:\Armando\LivemRNA\mRNADynamics\lib\audio\embryo_recorded_as_unhealthy.m4a';
    if exist(audioFile, 'file')
        [y, Fs] = audioread(audioFile);
    end
    
end

NFrames = liveExperiment.nFrames;
NSlices = liveExperiment.zDim;
yDim = liveExperiment.yDim; 
xDim = liveExperiment.xDim;

NChannels = size(movieMat, 5);


% initializes cell arrays for all the histone projections
median_proj = cell(NChannels, ceil(sum(NFrames) / skip_factor));
max_proj = cell(NChannels, ceil(sum(NFrames) / skip_factor));
middle_proj = cell(NChannels, ceil(sum(NFrames) / skip_factor));
midsum_proj = cell(NChannels, ceil(sum(NFrames) / skip_factor));
custom_proj = cell(NChannels, ceil(sum(NFrames) / skip_factor));


truncateAtColon = @(str) str(1:strfind(str, ':')-1);


if iscell(Channel1)
    Channel1 = Channel1{1};
end
if iscell(Channel2)
    Channel2 = Channel2{1};
end
if iscell(Channel3)
    Channel3 = Channel3{1};
end
ch1pre =  truncateAtColon(Channel1);
ch2pre =  truncateAtColon(Channel2);
ch3pre =  truncateAtColon(Channel3);

if contains(Channel1, ':')
    Channel1 = truncateAtColon(Channel1);
end
if contains(Channel2, ':')
    Channel2 = truncateAtColon(Channel2);
end
if contains(Channel3, ':')
    Channel3 = truncateAtColon(Channel3);
end

projectionChannels = {Channel1, Channel2, Channel3};

%construct cell to store projections for each frame separately
projCell = cell(NFrames, 1);
chCell = cell(NFrames, 1);
for ff = 1:NFrames
    projCell{ff} = ProjectionType;
    chCell{ff} = {Channel1, Channel2, Channel3};
end


% creates and stores histone slices
% idx = 1;
for framesIndex = 1:NFrames
    %         if mod(idx, skip_factor) == 1
    for channelIndex = 1:NChannels
        
        HisSlices = movieMat(:, :, :, framesIndex,channelIndex); %ch z t y x
        %
        %                 median_proj{channelIndex, ceil(idx / skip_factor)} = calculateProjection(...
        %                     'medianprojection', NSlices, HisSlices);
        %                 max_proj{channelIndex, ceil(idx / skip_factor)} = calculateProjection(...
        %                     'maxprojection', NSlices, HisSlices);
        %                 middle_proj{channelIndex, ceil(idx / skip_factor)} = calculateProjection(...
        %                     'middleprojection', NSlices, HisSlices);
        %                 midsum_proj{channelIndex, ceil(idx / skip_factor)} = calculateProjection(...
        %                     'midsumprojection', NSlices, HisSlices);
        %                 custom_proj{channelIndex, ceil(idx / skip_factor)} = calculateProjection(...
        %                     'customprojection', NSlices, HisSlices, max_custom, min_custom);
        
        %                  median_proj{channelIndex, framesIndex} = calculateProjection(...
        %                     'medianprojection', NSlices, HisSlices);
        median_proj{channelIndex, framesIndex} = calculateProjection(...
            'midsumprojection', NSlices, HisSlices);
        max_proj{channelIndex, framesIndex} = calculateProjection(...
            'maxprojection', NSlices, HisSlices);
        middle_proj{channelIndex, framesIndex} = calculateProjection(...
            'middleprojection', NSlices, HisSlices);
        midsum_proj{channelIndex, framesIndex} = calculateProjection(...
            'midsumprojection', NSlices, HisSlices);
%         custom_proj{channelIndex, framesIndex} = calculateProjection(...
%             'customprojection', NSlices, HisSlices, max_custom, min_custom);
    end
    
    %         end
    %      idx = idx + 1;
end

% numFrames = idx - 1;
numFrames = NFrames;
frame = 1;


% sets up channel dropdown options
options = cell(1, NChannels);
for j = 1:NChannels
    options{j} = ['Channel ' num2str(j)];
end

% lays out user interface
screen_size = get(0, 'screensize');
dim = [screen_size(3) * 0.6, screen_size(4) * 0.75];
dimVec = [dim(1), dim(2), dim(1), dim(2)]; %to easily normalize units
fig = uifigure('Position', [100, 100, dim(1), dim(2)], 'Name', 'Choose Histone Channels');
set(fig,'KeyPressFcn',@keycall)
imgAxis = uiaxes(fig, 'Position', [20, 20, dim(1) - 20, dim(2) * 0.5]);
hisPrecision = 'uint8';
blank = zeros(yDim, xDim, hisPrecision);
himage = imshow(blank, [], 'Parent', imgAxis);

frame_slider = uislider(fig, 'Limits', [1, numFrames], 'Value', 1, ...
    'Position', [dim(1) * 0.25, dim(2) * 0.6, dim(1) * 0.5, dim(2) * 0.1]);

frame_label = uilabel(fig,  'Text', ['Frame: ', num2str(frame)] , 'Position', ...
    [dim(1) * 0.5, dim(2) * 0.65, dim(1) * 0.1, dim(2) * 0.05]);

frame_slider.ValueChangedFcn = @updateHisImage;

%% display contrast stuff
maxPos = dimVec .* [.7, .2, .1, .1];
maxLabelPos = dimVec .* [.7, .225, .1, .05];
minPos = dimVec .* [.7, .4, .1, .1];
minLabelPos = dimVec .* [.7, .425, .1, .05];


max_label = uilabel(fig,  'Text', 'max display value', 'Position', ...
    maxLabelPos);
max_slider = uislider(fig, 'Limits', [0, 255], 'Value', 255, ...
    'Position', maxPos);
max_slider.ValueChangedFcn = @updateHisImage;


min_label = uilabel(fig,  'Text', 'min display value', 'Position', ...
    minLabelPos);
min_slider = uislider(fig, 'Limits', [0, 255], 'Value', 0, ...
    'Position',minPos);
min_slider.ValueChangedFcn = @updateHisImage;

%%
channel_label = uilabel(fig, 'Position', [10, dim(2) * 0.93, dim(1) * 0.125, dim(2) * 0.05], ...
    'Text', 'Channels');
channel_list = uilistbox(fig, 'Position', [10, dim(2) * 0.77, dim(1) * 0.125, dim(2) * 0.15], ...
    'MultiSelect', 'on', 'Items', options, ...
    'Value', {'Channel 1'});

channel_list.ValueChangedFcn = @updateHisImage;

invert_label = uilabel(fig, 'Position', [dim(1) * 0.2, dim(2) * 0.93, dim(1) * 0.125, dim(2) * 0.05], ...
    'Text', 'Inversions');
invert_list = uilistbox(fig, 'Position', [dim(1) * 0.2, dim(2) * 0.77, dim(1) * 0.125, dim(2) * 0.15], ...
    'MultiSelect', 'on', 'Items', options, ...
    'Value', {});

psa_label = uilabel(fig, 'Position', [10, dim(2) * .7, dim(1) * .5, dim(2) * .07], ...
    'Text', sprintf('%s \n %s', 'You can unselect a channel or','select multiple channels by holding down Ctrl'));

invert_list.ValueChangedFcn = @updateHisImage;

proj_type_label = uilabel(fig, 'Position', [dim(1) * 0.4, dim(2) * 0.93, dim(1) * 0.125, dim(2) * 0.05], ...
    'Text', 'Projection Type');
proj_type_dropdown = uidropdown(fig, 'Position', ...
    [dim(1) * 0.4, dim(2) * 0.85, dim(1) * 0.125, dim(2) * 0.08], ...
    'Items', {'maxprojection', 'medianprojection', 'middleprojection', 'midsumprojection', 'customprojection'}, ...
    'Value', {'midsumprojection'});

proj_type_dropdown.ValueChangedFcn = @updateHisImage;



%%prefix and folder saving area
if ~isempty(Prefix)
    
    saveFolderLabel = uilabel(fig, 'Position', [dim(1) * 0.3, dim(2) * 0.06, dim(1) * 0.5, dim(2) * 0.05], ...
        'Text', ['Project directory: ', PreProcPath, filesep]);
    %     saveFolderLabel = uilabel(fig, 'Position', [dim(1) * 0.3, dim(2) * 0.1, dim(1) * 0.5, dim(2) * 0.05], ...
    %                 'Text', ['Save folder: ',hisFolder]);
    %              sprintf('%s\n%s','Save folder: ',PreProcPath, filesep));
    
    PrefixLabel = uilabel(fig, 'Position', [dim(1) * 0.3, dim(2) * 0.01, dim(1) * 0.5, dim(2) * 0.05], ...
        'Text', ['Project: ', Prefix]);
    
end

% record anaphase frames


if ~isempty(Prefix)
    
    anaphaseTable = uitable(fig, 'Position', [dim(1) * 0.05, dim(2) * 0.2, dim(1) * 0.17, dim(2) * 0.27],...
        ...
        'ColumnName', [], 'ColumnEditable', true, 'DisplayDataChangedFcn', @tableUpdated,...
        ...
        'ColumnFormat', {'char', 'numeric'},...
        ...
        'RowName',{'nc9';'nc10'; 'nc11'; 'nc12';'nc13'; 'nc14'},...
        ...
        'Data',anaphaseFramesInitial);
    
    saveAnaphasesButton = uibutton(fig, 'Text', 'Save anaphase frames', 'Position', ...
        [dim(1) * 0.05, dim(2) * 0.08, dim(1) * 0.2, dim(2) * 0.05],...
        'ButtonPushedFcn', @saveAnaphasesButtonPushed);
    
end

%keyboard

keyboardButton = uibutton(fig, 'Text', 'keyboard', 'Position', ...
    [dim(1) * 0.85, dim(2) * .05, dim(1) * 0.1, dim(2) * 0.05],...
    'ButtonPushedFcn', @keyboardButtonPushed);


% Create a check box:
% cbx = uicheckbox(fig,'Position', [dim(1) * 0.85, dim(2) * .3, dim(1) * 0.2, dim(2) * 0.05],...
%     'ValueChangedFcn',@cBoxChanged, 'Text',...
%     sprintf('%s \n %s', 'Flag embryo', 'as unhealthy'));
cbx = uicheckbox(fig,'Position', [dim(1) * 0.85, dim(2) * .3, dim(1) * 0.2, dim(2) * 0.05],...
    'Text',...
    sprintf('%s \n %s', 'Flag embryo', 'as unhealthy'));

embryoHealthButton = uibutton(fig, 'Text', sprintf('%s \n %s', 'Save embryo', 'health status'),...
    'Position', [dim(1) * 0.85, dim(2) * .2, dim(1) * 0.15, dim(2) * 0.1],...
    'ButtonPushedFcn', @saveEmbryoHealthButtonPushed);




% custom_label = uilabel(fig, 'Position', [dim(1) * 0.4, dim(2) * 0.75, dim(1) * 0.125, dim(2) * 0.05],...
%     'Text', 'Averaging Range');
% custom_edit_text = uieditfield(fig, 'Position', [dim(1) * 0.525, dim(2) * 0.75, dim(1) * 0.125, dim(2) * 0.05],...
%     'Value', [num2str(max_custom) ':' num2str(min_custom)]);
%
% custom_change_confirmation = uibutton(fig, 'Text', 'Update Custom Projection', ...
%     'Position', [dim(1) * .7, dim(2) * 0.75, dim(1) * 0.2, dim(2) * 0.05]);
% custom_change_confirmation.ButtonPushedFcn = @updatedCustom;

save_button = uibutton(fig, 'Text', 'Save Channel Selection', 'Position', ...
    [dim(1) * 0.6, dim(2) * 0.85, dim(1) * 0.2, dim(2) * 0.08]);

save_button.ButtonPushedFcn = @saveOptions;

updateHisImage();
uiwait(fig);

% called whenever the frame, channels, or projection is changed. Updates
% the histone image the UI shows
    function updateHisImage(~, ~)
        channels_to_use = channel_list.Value;
        inverted_channels = invert_list.Value;
        projection_type = proj_type_dropdown.Value;
        frame_slider.Value = ((ceil(frame_slider.Value / skip_factor) - 1) * skip_factor) + 1;
        
        frame = ceil(frame_slider.Value / skip_factor);
        frame_label.Text=  ['Frame: ', num2str(frame)];
        
        channels = [];
        for i = 1:3
            if any(strcmp(channels_to_use, ['Channel ' num2str(i)]))
                channels = [channels i];
            end
        end
        ProjectionTemp = [];
        for i = 1:length(channels)
            cIndex = channels(i);
            if strcmpi(projection_type, 'medianprojection')
                ProjectionTemp(:, :, i) = median_proj{cIndex, frame};
            elseif strcmpi(projection_type, 'middleprojection')
                ProjectionTemp(:, :, i) = middle_proj{cIndex, frame};
            elseif strcmpi(projection_type, 'midsumprojection')
                ProjectionTemp(:,:, i) = midsum_proj{cIndex,frame};
            elseif strcmpi(projection_type, 'maxprojection')
                ProjectionTemp(:, :, i) = max_proj{cIndex, frame};
            else
                ProjectionTemp(:, :, i) = custom_proj{cIndex, frame};
            end
            if any(strcmp(inverted_channels, ['Channel ' num2str(cIndex)]))
                ProjectionTemp(:, :, i) = imcomplement(ProjectionTemp(:, :, i));
            end
            % Use the reference histogram to scale the Projection (This part
            % might need some more optimization later-YJK)
            ProjectionTemp(:, :, i) = histeq(mat2gray(ProjectionTemp(:, :, i)), ReferenceHist);
            ProjectionTemp(:, :, i) = ProjectionTemp(:, :, i) *255;
            
            
        end
        
        % Getting average of all Projections
        if length(channels) > 1
            Projection = nanmean(ProjectionTemp, 3);
        else
            Projection = ProjectionTemp;
        end
        
        maxDisplayIntensity = round(max_slider.Value);
        minDisplayIntensity = round(min_slider.Value);
        
        if minDisplayIntensity >= maxDisplayIntensity
            maxDisplayIntensity = max(max(Projection));
            minDisplayIntensity = median(median(Projection));
        end
        
        himage.CData = Projection;
        try imgAxis.CLim = [minDisplayIntensity, maxDisplayIntensity]; catch; end
        
        projCell{frame} = projection_type;
        %         chCell{frame} = getChannels;
        projectionChannels = retrieveChannels;
        Channel1 = projectionChannels{1}; 
        Channel2 = projectionChannels{2};
        Channel3 = projectionChannels{3};
        
        
    end

% closes UI and returns chosen options
    function saveOptions(~, ~)
        ProjectionType = proj_type_dropdown.Value;
        if strcmpi(ProjectionType, 'customprojection')
            ProjectionType = [ProjectionType ':' num2str(max_custom) ':' num2str(min_custom)];
        end
        
        projectionChannels = retrieveChannels;
        Channel1 = projectionChannels{1};
        Channel2 = projectionChannels{2}; 
        Channel3 = projectionChannels{3};
        
        if returnHisMat
            
            hisMat = zeros(yDim, xDim, sum(NFrames), hisPrecision); % y x f
            
            for f = 1:NFrames
                hisMat(:, :, f) = generateNuclearChannel2(ProjectionType,...
                    projectionChannels, ReferenceHist, movieMat, f);
            end
            
            saveNuclearProjection(hisMat, [liveExperiment.preFolder, filesep, Prefix, '-His.tif']);
            
        end
        
        
        
        save(projectionTypeFile,'ProjectionType','-v6')
        save(channelsFile,'projectionChannels','-v6')
        
        save(anaphaseFile, 'anaphaseFrames', '-v6')
        
        isUnhealthy = cbx.Value;
        save(isUnhealthyFile, 'isUnhealthy','-v6');
        
        close(fig);
        
    end

    function updatedCustom(~, ~)
        
        % updates hisslice range and forces range to be valid
        try
            custom_value = strsplit(custom_edit_text.Value, ':');
            max_custom = max(min(round(str2double(custom_value(1))), NSlices(1)), 1);
            min_custom = max(min(round(str2double(custom_value(2)), NSlices(1))), max_custom);
            custom_edit_text.Value = [num2str(max_custom) ':' num2str(min_custom)];
        catch
            min_custom = max(max_custom, min_custom);
            custom_edit_text.Value = [num2str(max_custom) ':' num2str(min_custom)];
        end
        
        %redoes projection for custom_proj
        idx2 = 1;
        for frameIndex = 1:NFrames
            if mod(idx2, skip_factor) == 1
                for chanIndex = 1:NChannels
                    
                    HisSlices = squeeze(movieMat(:, :, :, frameIndex,chanIndex)); %ch z t y x
                    
                    custom_proj{chanIndex, ceil(idx2 / skip_factor)} = calculateProjection(...
                        'customprojection', NSlices, HisSlices, max_custom, min_custom);
                end
            end
            idx2 = idx2 + 1;
        end
        
        updateHisImage();
    end

    function Channels = retrieveChannels()

        if contains(Channel1, ':')
            Channel1 = truncateAtColon(Channel1);
        end
        
        if any(strcmp(channel_list.Value, 'Channel 1'))
            if any(strcmp(invert_list.Value, 'Channel 1'))
                Channel1 = [ch1pre, ':invertedNuclear'];
            else
                Channel1 = [ch1pre, ':Nuclear'];
            end
        else
            Channel1 = '';
        end
        
        
        
        if contains(Channel2, ':')
            Channel2 = truncateAtColon(Channel2);
        end
        
        if any(strcmp(channel_list.Value, 'Channel 2'))
            if any(strcmp(invert_list.Value, 'Channel 2'))
                Channel2 = [ch2pre, ':invertedNuclear'];
            else
                Channel2 = [ch2pre, ':Nuclear'];
            end
        else
            Channel2 = '';
        end
        
        
        
        if contains(Channel3, ':')
            Channel3 = truncateAtColon(Channel3);
        end
        
        if any(strcmp(channel_list.Value, 'Channel 3'))
            if any(strcmp(invert_list.Value, 'Channel 3'))
                Channel3 = [ch3pre, ':invertedNuclear'];
            else
                Channel3 = [ch3pre, ':Nuclear'];
            end
        else
            Channel3 = '';
        end
        
        
        Channels = {Channel1, Channel2, Channel3};
        
    end




    function tableUpdated(src,event)
        anaphaseFrames = anaphaseTable.Data;
    end

    function keyboardButtonPushed(src,event)
        keyboard;
    end

    function saveAnaphasesButtonPushed(src,event)
        save(anaphaseFile, 'anaphaseFrames', '-v6')
        disp(num2str(anaphaseFrames'));
        disp('Anaphase frames saved.');
    end

    function saveEmbryoHealthButtonPushed(src,event)
        
        isUnhealthy = cbx.Value;
        save(isUnhealthyFile, 'isUnhealthy','-v6');
        
        if isUnhealthy
            disp('Embryo recorded as unhealthy.');
            nBits = 16;
            sound(y, Fs, nBits);
        else
            disp('Embryo recorded as healthy.');
        end
    end


    function keycall(h,e)
        
        if strcmpi(e.Key, 'rightarrow')
            if frame_slider.Value + 1 <= NFrames
                frame_slider.Value = frame_slider.Value + 1;
                updateHisImage;
            end
        elseif strcmpi(e.Key, 'leftarrow')
            if frame_slider.Value - 1 >= 1
                frame_slider.Value = frame_slider.Value - 1;
                updateHisImage;
            end
        end
        
    end


end