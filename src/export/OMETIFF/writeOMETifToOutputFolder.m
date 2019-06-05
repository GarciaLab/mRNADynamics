function writeOMETifToOutputFolder(TIFImages, Prefix, zSize, frameSize, ProjectionType, Channel1, Channel2, OutputFolder)

  % TO-DO refactor projection logic, it should be the same regardless of
  % the microscope. We can do this after we've finished supporting
  % ome-tiff. JP 2019-05-29.
  % Check how many channels have ":Nuclear" in the MovieDatabase.csv
  NuclearChannels = [contains(Channel1, 'Nuclear', 'IgnoreCase', true), contains(Channel2, 'Nuclear', 'IgnoreCase', true)];
  nNuclearChannels = sum(NuclearChannels);
  NuclearChannel = find(~ NuclearChannels == 0, nNuclearChannels);
    
  FramesFormat = '%03d';
  ChannelsFormat = '%02d';
  
  % temporary, this produces only one trailing zero,
  % so will produce '123' for '123' but '05' for '5', instead of '005'
  % we still have to check how to change z-stacks for more than 99 slices
  % without affecting existing code.
  ZSlicesFormat = '%02d'; 
  
  % we initialize a BlankImage to write at the start and end of each z-stack
  % with the same size as the first original image
  blankImageSize = size(TIFImages{1}{1,1});
  fprintf('BlankImage size will be %s\n', num2str(blankImageSize));
  BlankImage = uint16(zeros(blankImageSize));
  
  NPlanes = size(TIFImages{1},1);

  % generate tif files for non-nuclear channels
  for i = 1:NPlanes
    imageData = TIFImages{1}{i,1};
    
    imageMetaData = TIFImages{1}{i,2};
    imageMetaData = strtrim(split(imageMetaData, ';'));
    
    CurrentZSlice = str2double(getValueFromOMEMetadata(imageMetaData, 3));
    CurrentNChannel = getValueFromOMEMetadata(imageMetaData, 4);
    CurrentNFrame = getValueFromOMEMetadata(imageMetaData, 5);
    
    if str2double(CurrentNChannel) == NuclearChannel
      fprintf('Frame %s, Channel %d is nuclear, skipping tif generation. Will generate histone projection tif file instead.\n', CurrentNFrame, NuclearChannel);
      % add plane to NuclearPlane matrix to post-process later to generate
      % Histone projections
      NuclearPlanes{CurrentZSlice, str2double(CurrentNFrame)} = TIFImages{1}{i};
      continue;
    end
    
    if CurrentZSlice == 1
      % Save a blank image at the beginning of the z-stack
      FileName = buildOMEFileName(Prefix, FramesFormat, CurrentNFrame, ChannelsFormat, CurrentNChannel, ZSlicesFormat, num2str(CurrentZSlice));
    
      fprintf('Writing BlankImage at the beginning of the z-stack for frame %s, channel %s, zSlice %d to %s\n',...
        CurrentNFrame, CurrentNChannel, CurrentZSlice, [OutputFolder,filesep,FileName]);
    
      imwrite(BlankImage, [OutputFolder,filesep,FileName]);
    end
        
    % due to we adding blank images, the LivemRNA zSlice number is
    % off-by-one with the original sZlices numbers.
    CurrentZSlice = CurrentZSlice + 1;
        
    FileName = buildOMEFileName(Prefix, FramesFormat, CurrentNFrame, ChannelsFormat, CurrentNChannel, ZSlicesFormat, num2str(CurrentZSlice));
    
    fprintf('Writing TIF file for frame %s, channel %s, zSlice %d to %s\n', CurrentNFrame, CurrentNChannel, CurrentZSlice, [OutputFolder,filesep,FileName]);
    imwrite(imageData, [OutputFolder,filesep,FileName]);

    if CurrentZSlice > zSize
      % Save a blank image at the end of the z-stack
      CurrentZSlice = CurrentZSlice + 1;
      FileName = buildOMEFileName(Prefix, FramesFormat, CurrentNFrame, ChannelsFormat, CurrentNChannel, ZSlicesFormat, num2str(CurrentZSlice));
    
      fprintf('Writing BlankImage at the end of the z-stack for frame %s, channel %s, zSlice %d to %s\n',...
        CurrentNFrame, CurrentNChannel, CurrentZSlice, [OutputFolder,filesep,FileName]);
    
      imwrite(BlankImage, [OutputFolder,filesep,FileName]);
    end
  end
  
  
  % calculate z-stack projection for each frame
  Projections = calculateOMETIFFProjections(ProjectionType, NuclearPlanes);
  
  % write projection files
  for i = 1:frameSize
    imwrite(uint16(Projections{i}), [OutputFolder, filesep, Prefix, '-His_', sprintf(FramesFormat, i), '.tif']);
  end
  
end
