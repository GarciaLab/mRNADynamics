% Edited from processFiducialChannel.m
% by Yang Joon Kim (yjkim90@berkeley.edu) on 10/23/2018
% This will check all channels, then grab channels with ":Nuclear" in the
% string, then add the Projected images from multiple channels.
% For MovieDatabase Channels, we should put ":Nuclear" or
% "invertedNuclear" for those channels to be recognized for histone
% channel generation.
function Projection = generateNuclearChannel(numberOfFrames, LIFImages, framesIndex, seriesIndex, NSlices, NChannels, nuclearChannel, ProjectionType, ExperimentType, Channel1, Channel2, Channel3, ReferenceHist, OutputFolder, Prefix)

  % Check how many channels have ":Nuclear" in the MovieDatabase.csv
  NuclearChannels = [contains(Channel1, 'Nuclear', 'IgnoreCase', true), ...
                       contains(Channel2, 'Nuclear', 'IgnoreCase', true), ...
                       contains(Channel3, 'Nuclear', 'IgnoreCase', true)];
  nNuclearChannels = sum(NuclearChannels);

  InvertedChannels = [contains(Channel1, 'inverted', 'IgnoreCase', true), ...
                        contains(Channel2, 'inverted', 'IgnoreCase', true), ...
                        contains(Channel3, 'inverted', 'IgnoreCase', true)];
  
  if nNuclearChannels ~= 0

    for ChannelIndex = 1:nNuclearChannels

      % Find the corresponding channel
      AllNuclearChannels = find(~ NuclearChannels == 0, nNuclearChannels);
      nuclearChannel = AllNuclearChannels(ChannelIndex);

      % For all 'nuclear' channels, generate HisSlices, and do projection
      HisSlices = generateHisSlices(LIFImages, NSlices, NChannels, nuclearChannel, framesIndex, seriesIndex);

      ProjectionTemp(:, :, ChannelIndex) = calculateProjection(ProjectionType, NSlices, HisSlices);

      % Think about "invertedNuclear", for example, MCP-mCherry, then
      % invert the ProjectionTemp using imcomplement
      if InvertedChannels(nuclearChannel) == 1
        ProjectionTemp(:, :, ChannelIndex) = imcomplement(ProjectionTemp(:, :, ChannelIndex));
      end

      % Use the reference histogram to scale the Projection (This part
      % might need some more optimization later-YJK)
      ProjectionTemp(:, :, ChannelIndex) = histeq(mat2gray(ProjectionTemp(:, :, ChannelIndex)), ReferenceHist);
      ProjectionTemp(:, :, ChannelIndex) = ProjectionTemp(:, :, ChannelIndex) * 10000;
    end

    % Getting average of all Projections
    if nNuclearChannels > 1
      Projection = nanmean(ProjectionTemp, 3);
    elseif nNuclearChannels == 1
      Projection = ProjectionTemp;
    end

    % In case of old datasets (does not have ":Nuclear")
  else

    HisSlices = generateHisSlices(LIFImages, NSlices, NChannels, nuclearChannel, framesIndex, seriesIndex);

    Projection = calculateProjection(ProjectionType, NSlices, HisSlices);
    
    %YJK : Think about the case when there is no His channel,
    %and it is inputoutput mode or 1spot mode or 2spot2color.
    %We can use (MCP-mCherry) either inverted or raw
    %images to make fake histone images.
    if (isempty(strfind(Channel1{1}, 'His'))) && (isempty(strfind(Channel2{1}, 'His'))) && (isempty(strfind(Channel3{1}, 'His')))
      
      if strcmpi(ExperimentType, 'inputoutput') | strcmpi(ExperimentType, '1spot') | strcmpi(ExperimentType, '2spot2color') | strcmpi(ExperimentType, 'input')
        
        if (~ isempty(strfind(Channel1{1}, 'NLS'))) | (~ isempty(strfind(Channel2{1}, 'NLS')))
          %don't invert with NLS-MCP-mCherry
        else
          %We don't want to use all slices. Only the center ones
          StackCenter = round((min(NSlices) - 1) / 2);
          StackRange = StackCenter - 1:StackCenter + 1;
          
          if strcmp(ProjectionType, 'medianprojection')
            Projection = median(HisSlices(:, :, StackRange), [], 3);
          else
            Projection = max(HisSlices(:, :, StackRange), [], 3);
          end
          
          %invert images to make nuclei bright
          Projection = imcomplement(Projection);
        end
        
        Projection = histeq(mat2gray(Projection), ReferenceHist);
        Projection = Projection * 10000;
      end
      
    end
    
  end
  
  imwrite(uint16(Projection), [OutputFolder, filesep, Prefix, '-His_', iIndex(numberOfFrames, 3), '.tif']);
  
end

function HisSlices = generateHisSlices(LIFImages, NSlices, NChannels, fiducialChannel, framesIndex, seriesIndex)
  
  % For all 'nuclear' channels, generate HisSlices, and do projection
  HisSlices = zeros([size(LIFImages{seriesIndex}{1, 1}, 1), size(LIFImages{seriesIndex}{1, 1}, 2), NSlices(seriesIndex)]);
  z = 1;
  firstImage = (framesIndex - 1) * NSlices(seriesIndex) * NChannels + 1 + (fiducialChannel - 1);
  lastImage = framesIndex * NSlices(seriesIndex) * NChannels;
  
  for imagesIndex = firstImage:NChannels:lastImage
    HisSlices(:, :, z) = LIFImages{seriesIndex}{imagesIndex, 1};
    z = z + 1;
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
