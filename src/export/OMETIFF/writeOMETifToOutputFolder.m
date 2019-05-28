function writeOMETifToOutputFolder(TIFImages, Prefix, OutputFolder)
  FramesFormat = '%03d';
  ChannelsFormat = '%02d';
  
  % temporary, this produces only one trailing zero,
  % so will produce '123' for '123' but '05' for '5', instead of '005'
  % we still have to check how to change z-stacks for more than 99 slices
  % without affecting existing code.
  ZSlicesFormat = '%02d'; 

  NPlanes = size(TIFImages{1},1);

  for i = 1:NPlanes
    imageData = TIFImages{1}{i,1};
    
    imageMetaData = TIFImages{1}{i,2};
    imageMetaData = strtrim(split(imageMetaData, ';'));
    
    CurrentZSlice = getValueFromOMEMetadata(imageMetaData, 3);
    
    
    CurrentNChannel = getValueFromOMEMetadata(imageMetaData, 4);
    CurrentNFrame = getValueFromOMEMetadata(imageMetaData, 5);
    
    
    FileName = buildOMEFileName(Prefix, FramesFormat, CurrentNFrame, ChannelsFormat, CurrentNChannel, ZSlicesFormat, CurrentZSlice);
    fprintf('Writing TIF file for frame %s, channel %s, zSlice %s to %s\n', CurrentNFrame, CurrentNChannel, CurrentZSlice, [OutputFolder,filesep,FileName]);
    
    imwrite(imageData, [OutputFolder,filesep,FileName]);
  end
end
