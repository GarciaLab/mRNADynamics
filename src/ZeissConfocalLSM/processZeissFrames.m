function processZeissFrames(Prefix, ExperimentType, Channel1, Channel2, OutputFolder, LSMImages, LSMIndex, FrameRange, NSlices, NChannels, coatChannel, fiducialChannel, ReferenceHist)
  % Create a blank image
  BlankImage = uint16(zeros(size(LSMImages{1}{1, 1})));

  % Redefine LSMImages as a cell. I think we need this so that the
  % parfor loop doesn't freak out.
  LSMImages = LSMImages{1};

  for framesIndex = 1:length(FrameRange)
    % First do the coat protein channel
    % Save the blank images at the beginning and end of the stack
    NameSuffix = ['_ch', iIndex(coatChannel, 2)];

    NewName = [Prefix, '_', iIndex(FrameRange(framesIndex), 3), '_z', iIndex(1, 2), NameSuffix, '.tif'];
    imwrite(BlankImage, [OutputFolder, filesep, NewName]);
    NewName = [Prefix, '_', iIndex(FrameRange(framesIndex), 3), '_z', iIndex(min(NSlices) + 2, 2), NameSuffix, '.tif'];
    imwrite(BlankImage, [OutputFolder, filesep, NewName]);
    
    %Copy the rest of the images
    slicesCounter = 1; %Counter for slices
    
    for slicesIndex = ((framesIndex - 1) * NSlices(LSMIndex) * NChannels(LSMIndex) + 1 + (coatChannel - 1)):...
      NChannels:(framesIndex * NSlices(LSMIndex)) * NChannels
      if slicesCounter <= min(NSlices)
        NewName = [Prefix, '_', iIndex(FrameRange(framesIndex), 3), '_z', iIndex(slicesCounter + 1, 2), NameSuffix, '.tif'];
        imwrite(LSMImages{slicesIndex, 1}, [OutputFolder, filesep, NewName]);
        slicesCounter = slicesCounter + 1;
      end
    end

    %Now do His-RFP
    HisSlices = zeros([size(LSMImages{1, 1}, 1), size(LSMImages{1, 1}, 2), NSlices(LSMIndex)]);
    slicesCounter = 1;
    
    for slicesIndex = ((framesIndex - 1) * NSlices(LSMIndex) * NChannels(LSMIndex) + 1 + (fiducialChannel - 1)):NChannels(LSMIndex):...
      (framesIndex * NSlices(LSMIndex)) * NChannels(LSMIndex)
      
      HisSlices(:, :, slicesCounter) = LSMImages{slicesIndex, 1};
      slicesCounter = slicesCounter + 1;
    end

    Projection = median(HisSlices, 3);
    if strcmpi(ExperimentType, 'inputoutput')
      %YJK : Think about the case when there is no His channel,
      %and it is inputoutput mode or 1spot mode or 2spot2color.
      %We can use (MCP-mCherry) either inverted or raw
      %images to make fake histone images.
      if (isempty(strfind(Channel1{1}, 'His'))) && (isempty(strfind(Channel2{1}, 'His')))
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

    imwrite(uint16(Projection), [OutputFolder, filesep, Prefix, '-His_', iIndex(FrameRange(framesIndex), 3), '.tif']);
  end

end
