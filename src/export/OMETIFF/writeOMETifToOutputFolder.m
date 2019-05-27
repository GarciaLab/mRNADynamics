function writeOMETifToOutputFolder(TIFImages, Prefix, NFrames, NChannels, NZSlices, OutputFolder)
  FramesFormat = getStringFormatForFileName(NFrames);
  ChannelsFormat = getStringFormatForFileName(NChannels);
  ZSlicesFormat = getStringFormatForFileName(NZSlices);

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
