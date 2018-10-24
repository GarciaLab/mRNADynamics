function Projection = generateNuclearChannel(numberOfFrames, LIFImages, framesIndex, seriesIndex, NSlices, NChannels, fiducialChannel, histoneChannel, ProjectionType, ExperimentType, Channel1, Channel2, Channel3, ReferenceHist, OutputFolder, Prefix)
  % Edited from processFiducialChannel.m 
  % by Yang Joon Kim (yjkim90@berkeley.edu) on 10/23/2018
  % This will check all channels, then grab channels with ":Nuclear" in the
  % string, then add the Projected images from multiple channels.
  % For MovieDatabase Channels, we should put ":Nuclear" or
  % "invertedNuclear" for those channels to be recognized for histone
  % channel generation.
  
  % Check how many channels have ":Nuclear" in the MovieDatabase.csv
  NuclearChannels = [contains(Channel1,'Nuclear','IgnoreCase',true),...
                contains(Channel2,'Nuclear','IgnoreCase',true),...
                contains(Channel3,'Nuclear','IgnoreCase',true)];
  nNuclearChannels = sum(NuclearChannels);

  InvertedChannels = [contains(Channel1,'inverted','IgnoreCase',true),...
                contains(Channel2,'inverted','IgnoreCase',true),...
                contains(Channel3,'inverted','IgnoreCase',true)];
  
                    
  if ~nNuclearChannels==0
      for nuclearChannel=1:nNuclearChannels
          
        % Find the corresponding channel
        AllNuclearChannels = find(~NuclearChannels==0,nNuclearChannels);
        nuclearChannel = AllNuclearChannels(nuclearChannel);
        
        % For all 'nuclear' channels, generate HisSlices, and do projection
        HisSlices = zeros([size(LIFImages{seriesIndex}{1,1},1), size(LIFImages{seriesIndex}{1,1},2), NSlices(seriesIndex)]);
        n = 1;
        firstImage = (framesIndex-1) * NSlices(seriesIndex) * NChannels + 1 + (nuclearChannel - 1);
        lastImage = framesIndex * NSlices(seriesIndex) * NChannels;
        for imagesIndex = firstImage:NChannels:lastImage
            HisSlices(:,:,n) = LIFImages{seriesIndex}{imagesIndex,1};
            n = n + 1;
        end
        
        % Calculate the projection (either Maximum or Median)
          if strcmpi(ProjectionType,'medianprojection')
            ProjectionTemp(:,:,nuclearChannel) = median(HisSlices, 3);
          elseif strcmpi(ProjectionType,'middleprojection')
            ProjectionTemp(:,:,nuclearChannel) = max(HisSlices(:,:,round(NSlices*.50):round(NSlices*.75)), [], 3); 
          else
            ProjectionTemp(:,:,nuclearChannel) = max(HisSlices, [], 3);
          end
        % Think about "invertedNuclear", for example, MCP-mCherry, then
        % invert the ProjectionTemp using imcomplement
        if InvertedChannels(nuclearChannel)==1
            ProjectionTemp(:,:,nuclearChannel) = imcomplement(ProjectionTemp(:,:,nuclearChannel));
        end
        
        % Use the reference histogram to scale the Projection (This part
        % might need some more optimization later-YJK)
        ProjectionTemp(:,:,nuclearChannel) = histeq(mat2gray(ProjectionTemp(:,:,nuclearChannel)), ReferenceHist);
        ProjectionTemp(:,:,nuclearChannel) = ProjectionTemp(:,:,nuclearChannel) * 10000;   
      end
    
      Projection = nanmean(ProjectionTemp,3);
  % In case of old datasets (does not have ":Nuclear")
  else
        HisSlices = zeros([size(LIFImages{seriesIndex}{1,1},1), size(LIFImages{seriesIndex}{1,1},2), NSlices(seriesIndex)]);
        n = 1;
        firstImage = (framesIndex-1) * NSlices(seriesIndex) * NChannels + 1 + (fiducialChannel - 1);
        lastImage = framesIndex * NSlices(seriesIndex) * NChannels;
        
        for imagesIndex = firstImage:NChannels:lastImage
            HisSlices(:,:,n) = LIFImages{seriesIndex}{imagesIndex,1};
            n = n + 1;
        end
      
      % Projection
      if strcmpi(ProjectionType,'medianprojection')
        Projection = median(HisSlices, 3);
      elseif strcmpi(ProjectionType,'middleprojection')
        Projection = max(HisSlices(:,:,round(NSlices*.50):round(NSlices*.75)), [], 3); 
      else
        Projection = max(HisSlices, [], 3);
      end

      %YJK : Think about the case when there is no His channel,
      %and it is inputoutput mode or 1spot mode or 2spot2color.
      %We can use (MCP-mCherry) either inverted or raw
      %images to make fake histone images.
      if (isempty(strfind(Channel1{1}, 'His'))) && (isempty(strfind(Channel2{1}, 'His'))) && (isempty(strfind(Channel3{1}, 'His')))
        if strcmpi(ExperimentType, 'inputoutput')|strcmpi(ExperimentType, '1spot')|strcmpi(ExperimentType,'2spot2color')|strcmpi(ExperimentType,'input')
          if (~isempty(strfind(Channel1{1}, 'NLS')))|(~isempty(strfind(Channel2{1}, 'NLS')))
            %don't invert with NLS-MCP-mCherry
          else
            %We don't want to use all slices. Only the center ones
            StackCenter = round((min(NSlices) - 1) / 2);
            StackRange = StackCenter - 1:StackCenter + 1;
            if strcmp(ProjectionType, 'medianprojection')
                Projection = median(HisSlices(:,:,StackRange), [], 3);
            else
                Projection = max(HisSlices(:,:,StackRange), [], 3);
            end
            %invert images to make nuclei bright
            Projection = imcomplement(Projection);
          end
          Projection = histeq(mat2gray(Projection), ReferenceHist);
          Projection = Projection * 10000;
        end
      end   
  end


  imwrite(uint16(Projection),[OutputFolder, filesep, Prefix, '-His_', iIndex(numberOfFrames, 3), '.tif']);
  
end
