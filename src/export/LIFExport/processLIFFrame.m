function processLIFFrame(numberOfFrames, Prefix, BlankImage, OutputFolder, LIFImages, framesIndex,...
  seriesIndex, NChannels, NSlices, ExperimentType, Channel1, Channel2, Channel3, ProjectionType,...
  fiducialChannel, histoneChannel, ReferenceHist, coatChannel, inputProteinChannel, zslicesPadding)
  for channelIndex = 1:NChannels
    processLIFChannel(ExperimentType, channelIndex, numberOfFrames, Prefix, BlankImage, OutputFolder,...
      LIFImages, framesIndex, seriesIndex, NChannels, NSlices, coatChannel, inputProteinChannel, zslicesPadding);
  end
  %Now copy nuclear tracking images
  if fiducialChannel
    generateNuclearChannel(numberOfFrames, LIFImages, framesIndex, seriesIndex, NSlices, NChannels, fiducialChannel, ProjectionType, ExperimentType, Channel1, Channel2, Channel3, ReferenceHist, OutputFolder, Prefix);
  end
end
