function processZeissFrames(Prefix, OutputFolder, LSMImages, LSMIndex, FrameRange, NSlices, NChannels, coatChannel, fiducialChannel)
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
    imwrite(uint16(Projection), [OutputFolder, filesep, Prefix, '-His_', iIndex(FrameRange(framesIndex), 3), '.tif']);
  end

end
