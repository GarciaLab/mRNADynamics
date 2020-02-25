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

function [Channels, ProjectionType] = chooseNuclearChannels(...
    LIFImages, NSeries, NSlices, NChannels, NFrames, ProjectionType, Channels, ReferenceHist)

skip_factor = 5; % Only uses 1/skip_factor frames

% initializes cell arrays for all the histone projections
median_proj = cell(NChannels, ceil(sum(NFrames) / skip_factor));
max_proj = cell(NChannels, ceil(sum(NFrames) / skip_factor));
middle_proj = cell(NChannels, ceil(sum(NFrames) / skip_factor));
custom_proj = cell(NChannels, ceil(sum(NFrames) / skip_factor));
mean_proj = cell(NChannels, ceil(sum(NFrames) / skip_factor));

Channel1 = Channels{1}; Channel2 = Channels{2}; Channel3 = Channels{3};

% default custom projection parameters
max_custom = 1; % highest histone channel slice used
min_custom = 5; % lowest histone channel slice used

% creates and stores histone slices
idx = 1;
for seriesIndex = 1:NSeries
    for framesIndex = 1:NFrames(seriesIndex)
        if mod(idx, skip_factor) == 1
            for channelIndex = 1:NChannels
                HisSlices = generateHisSlices(LIFImages, NSlices, NChannels, ...
                    channelIndex, framesIndex, seriesIndex);
                median_proj{channelIndex, ceil(idx / skip_factor)} = calculateProjection(...
                    'medianprojection', NSlices, HisSlices);
                max_proj{channelIndex, ceil(idx / skip_factor)} = calculateProjection(...
                    'maxprojection', NSlices, HisSlices);
                middle_proj{channelIndex, ceil(idx / skip_factor)} = calculateProjection(...
                    'middleprojection', NSlices, HisSlices);
                custom_proj{channelIndex, ceil(idx / skip_factor)} = calculateProjection(...
                    'customprojection', NSlices, HisSlices, max_custom, min_custom);
            end       
        end
        idx = idx + 1;
    end
end
numFrames = idx - 1;

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
img = uiaxes(fig, 'Position', [20, 20, dim(1) - 20, dim(2) * 0.5]);

frame_label = uilabel(fig,  'Text', 'Frame', 'Position', ...
    [dim(1) * 0.5, dim(2) * 0.65, dim(1) * 0.1, dim(2) * 0.05]);

frame_slider = uislider(fig, 'Limits', [1, numFrames], 'Value', 1, ...
    'Position', [dim(1) * 0.25, dim(2) * 0.6, dim(1) * 0.5, dim(2) * 0.1]);

frame_slider.ValueChangedFcn = @updateHisImage;

%% display contrast stuff
maxPos = dimVec .* [.7, .2, .1, .1];
maxLabelPos = dimVec .* [.7, .225, .1, .05];
minPos = dimVec .* [.7, .4, .1, .1];
minLabelPos = dimVec .* [.7, .425, .1, .05];

max_label = uilabel(fig,  'Text', 'max display value', 'Position', ...
    maxLabelPos);
max_slider = uislider(fig, 'Limits', [1, 10000], 'Value', 10000, ...
    'Position', maxPos);
max_slider.ValueChangedFcn = @updateHisImage;


min_label = uilabel(fig,  'Text', 'min display value', 'Position', ...
    minLabelPos);
min_slider = uislider(fig, 'Limits', [0, 10000], 'Value', 1, ...
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
    'Text', 'You can unselect a channel or select multiple channels by holding down Ctrl');

invert_list.ValueChangedFcn = @updateHisImage;

proj_type_label = uilabel(fig, 'Position', [dim(1) * 0.4, dim(2) * 0.93, dim(1) * 0.125, dim(2) * 0.05], ...
    'Text', 'Projection Type');
proj_type_dropdown = uidropdown(fig, 'Position', ...
    [dim(1) * 0.4, dim(2) * 0.85, dim(1) * 0.125, dim(2) * 0.08], ...
    'Items', {'maxprojection', 'medianprojection', 'middleprojection', 'customprojection'}, ...
    'Value', {'maxprojection'});

proj_type_dropdown.ValueChangedFcn = @updateHisImage;

custom_label = uilabel(fig, 'Position', [dim(1) * 0.4, dim(2) * 0.75, dim(1) * 0.125, dim(2) * 0.05],...
    'Text', 'Averaging Range');
custom_edit_text = uieditfield(fig, 'Position', [dim(1) * 0.525, dim(2) * 0.75, dim(1) * 0.125, dim(2) * 0.05],...
    'Value', [num2str(max_custom) ':' num2str(min_custom)]);

custom_change_confirmation = uibutton(fig, 'Text', 'Update Custom Projection', ...
    'Position', [dim(1) * .7, dim(2) * 0.75, dim(1) * 0.2, dim(2) * 0.05]);
custom_change_confirmation.ButtonPushedFcn = @updatedCustom;

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
            elseif strcmpi(projection_type, 'meanprojection')
                ProjectionTemp(:,:, i) = mean_proj(cIndex,frame);
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
              ProjectionTemp(:, :, i) = ProjectionTemp(:, :, i) * 10000;
        end

        % Getting average of all Projections
        if length(channels) > 1
          Projection = nanmean(ProjectionTemp, 3);
        else
          Projection = ProjectionTemp;
        end
          
        maxB = round(max_slider.Value);
        minB = round(min_slider.Value);
        
         if minB >= maxB
            maxB = max(max(Projection));
            minB = median(median(Projection));
        end
%         imshow(uint16(Projection), [median(median(Projection)), max(max(Projection))], 'Parent', img);
        imshow(uint16(Projection), [minB, maxB], 'Parent', img);
    end
    
    % closes UI and returns chosen options
    function saveOptions(~, ~)
        if ~isempty(Channel1)
            if any(strcmp(channel_list.Value, 'Channel 1'))
                if any(strcmp(invert_list.Value, 'Channel 1'))
                    Channel1{1} = [Channel1{1}, ':invertedNuclear'];
                else
                    Channel1{1} = [Channel1{1}, ':Nuclear'];
                end
            else
                Channel1 = '';
            end
        end
        if ~isempty(Channel2)
            if any(strcmp(channel_list.Value, 'Channel 2'))
                if any(strcmp(invert_list.Value, 'Channel 2'))
                    Channel2{1} = [Channel2{1}, ':invertedNuclear'];
                else
                    Channel2{1} = [Channel2{1}, ':Nuclear'];
                end
            else
                Channel2 = '';
            end
            
        end
        if ~isempty(Channel3)
            if any(strcmp(channel_list.Value, 'Channel 3'))
                if any(strcmp(invert_list.Value, 'Channel 3'))
                    Channel3{1} = [Channel3{1}, ':invertedNuclear'];
                else
                    Channel3{1} = [Channel3{1}, ':Nuclear'];
                end
            else
                Channel3 = '';
            end
        end
        
        Channels = {Channel1, Channel2, Channel3};
        ProjectionType = proj_type_dropdown.Value;
        if strcmpi(ProjectionType, 'customprojection')
            ProjectionType = [ProjectionType ':' num2str(max_custom) ':' num2str(min_custom)];
        end
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
        for seriesIndex = 1:NSeries
            for framesIndex = 1:NFrames(seriesIndex)
                if mod(idx2, skip_factor) == 1
                    for channelIndex = 1:NChannels
                        HisSlices = generateHisSlices(LIFImages, NSlices, NChannels, ...
                            channelIndex, framesIndex, seriesIndex);
                        custom_proj{channelIndex, ceil(idx2 / skip_factor)} = calculateProjection(...
                            'customprojection', NSlices, HisSlices, max_custom, min_custom);
                    end       
                end
                idx2 = idx2 + 1;
            end
        end
        updateHisImage();
    end

end