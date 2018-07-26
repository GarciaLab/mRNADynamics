function Projection = processFiducialChannel(numberOfFrames, LIFImages, framesIndex, seriesIndex, NSlices, NChannels, fiducialChannel, histoneChannel, ProjectionType, ExperimentType, Channel1, Channel2, Channel3, ReferenceHist, OutputFolder, Prefix)
  HisSlices = zeros([size(LIFImages{seriesIndex}{1,1},1), size(LIFImages{seriesIndex}{1,1},2), NSlices(seriesIndex)]);
  otherSlices = HisSlices;
  n = 1;
  firstImage = (framesIndex-1) * NSlices(seriesIndex) * NChannels + 1 + (fiducialChannel - 1);
  lastImage = framesIndex * NSlices(seriesIndex) * NChannels;
  for imagesIndex = firstImage:NChannels:lastImage
      HisSlices(:,:,n) = LIFImages{seriesIndex}{imagesIndex,1};
      otherSlices(:,:,n) = LIFImages{seriesIndex}{imagesIndex+1,1};
      n = n + 1;
  end

  if histoneChannel
    Projection = getHistoneChannelProjection(ProjectionType, HisSlices, ExperimentType, Channel1, Channel2, Channel3, NSlices, ReferenceHist);
  else 
    Projection = getDefaultChannelProjection(NSlices, HisSlices);
  end

  imwrite(uint16(Projection),[OutputFolder, filesep, Prefix, '-His_', iIndex(numberOfFrames, 3), '.tif']);
end
