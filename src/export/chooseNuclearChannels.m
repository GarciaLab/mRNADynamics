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

function [Channel1, Channel2, Channel3, ProjectionType] = chooseNuclearChannels(...
    LIFImages, NSeries, NSlices, NChannels, NFrames, ProjectionType, Channel1, Channel2, ...
    Channel3, ReferenceHist)

% generates all the HisImages
skip_factor = 4; % Only uses 1/skip_factor frames
median_proj = cell(NChannels, ceil(sum(NFrames) / skip_factor));
max_proj = cell(NChannels, sum(NFrames));
middle_proj = cell(NChannels, sum(NFrames));
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
fig = uifigure('Position', [100, 100, dim(1), dim(2)], 'Name', 'Choose Histone Channels');
img = uiaxes(fig, 'Position', [20, 20, dim(1) - 20, dim(2) * 0.5]);

frame_label = uilabel(fig,  'Text', 'Frame', 'Position', ...
    [dim(1) * 0.5, dim(2) * 0.65, dim(1) * 0.1, dim(2) * 0.05]);

frame_slider = uislider(fig, 'Limits', [1, numFrames], 'Value', 1, ...
    'Position', [dim(1) * 0.25, dim(2) * 0.6, dim(1) * 0.5, dim(2) * 0.1]);

frame_slider.ValueChangedFcn = @updateHisImage;

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
    'Items', {'maxprojection', 'medianprojection', 'middleprojection'}, ...
    'Value', {'maxprojection'});

proj_type_dropdown.ValueChangedFcn = @updateHisImage;

save_button = uibutton(fig, 'Text', 'Save Channel Selection', 'Position', ...
    [dim(1) * 0.6, dim(2) * 0.85, dim(1) * 0.2, dim(2) * 0.08]);

save_button.ButtonPushedFcn = @saveOptions;

updateHisImage();
uiwait(fig);

    function updateHisImage(~, ~)
        channels_to_use = channel_list.Value;
        inverted_channels = invert_list.Value;
        projection_type = proj_type_dropdown.Value;
        frame_slider.Value = ((ceil(frame_slider.Value / 4) - 1) * 4) + 1;
        frame = ceil(frame_slider.Value / 4);
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
            else
                ProjectionTemp(:, :, i) = max_proj{cIndex, frame};
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
        imshow(uint16(Projection), [], 'Parent', img);
    end

    function saveOptions(~, ~)
        if ~isempty(Channel1)
            if any(strcmp(channel_list.Value, 'Channel 1'))
                if any(strcmp(invert_list.Value, 'Channel 1'))
                    Channel1{1} = [Channel1{1}, ':invertedNuclear'];
                else
                    Channel1{1} = [Channel1{1}, ':Nuclear'];
                end
            end
        end
        if ~isempty(Channel2)
            if any(strcmp(channel_list.Value, 'Channel 2'))
                if any(strcmp(invert_list.Value, 'Channel 2'))
                    Channel2{1} = [Channel2{1}, ':invertedNuclear'];
                else
                    Channel2{1} = [Channel2{1}, ':Nuclear'];
                end
            end
        end
        if ~isempty(Channel3)
            if any(strcmp(channel_list.Value, 'Channel 3'))
                if any(strcmp(invert_list.Value, 'Channel 3'))
                    Channel3{1} = [Channel3{1}, ':invertedNuclear'];
                else
                    Channel3{1} = [Channel3{1}, ':Nuclear'];
                end
            end
        end
        
        ProjectionType = proj_type_dropdown.Value;
        close(fig);
    end

end

function HisSlices = generateHisSlices(LIFImages, NSlices, NChannels, fiducialChannel, framesIndex, seriesIndex)
  
  % For all 'nuclear' channels, generate HisSlices, and do projection
  HisSlices = zeros([size(LIFImages{seriesIndex}{1, 1}, 1), size(LIFImages{seriesIndex}{1, 1}, 2), NSlices(seriesIndex)]);
  n = 1;
  firstImage = (framesIndex - 1) * NSlices(seriesIndex) * NChannels + 1 + (fiducialChannel - 1);
  lastImage = framesIndex * NSlices(seriesIndex) * NChannels;
  
  for imagesIndex = firstImage:NChannels:lastImage
    HisSlices(:, :, n) = LIFImages{seriesIndex}{imagesIndex, 1};
    n = n + 1;
  end
  
end

function Projection = calculateProjection(ProjectionType, NSlices, HisSlices)
  % Calculate the projection (either Maximum or Median)
  if strcmpi(ProjectionType, 'medianprojection')
    Projection = median(HisSlices, 3);
  elseif strcmpi(ProjectionType, 'middleprojection')
    Projection = max(HisSlices(:, :, round(NSlices * .50):round(NSlices * .75)), [], 3);
  else
    Projection = max(HisSlices, [], 3);
  end

end


