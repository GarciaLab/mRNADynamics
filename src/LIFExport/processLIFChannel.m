function processLIFChannel(ExperimentType, channelIndex, numberOfFrames, Prefix, BlankImage, OutputFolder, LIFImages, framesIndex, seriesIndex, NChannels, NSlices, coatChannel, inputProteinChannel)
  experimentType1 = (strcmpi(ExperimentType,'1spot') || strcmp(ExperimentType,'2spot') || strcmp(ExperimentType,'2spot1color')) && channelIndex == coatChannel;
  experimentType2 = strcmpi(ExperimentType,'2spot2color') || strcmpi(ExperimentType, 'inputoutput');
  experimentType3 = strcmpi(ExperimentType, 'input') && sum(channelIndex == inputProteinChannel);

  %Save the blank images at the beginning and end of the stack
  NameSuffix = ['_ch',iIndex(channelIndex,2)];
  NewName = [Prefix, '_', iIndex(numberOfFrames,3), '_z', iIndex(1,2), NameSuffix, '.tif'];
  imwrite(BlankImage, [OutputFolder, filesep, NewName]);
  NewName = [Prefix, '_', iIndex(numberOfFrames,3), '_z', iIndex(min(NSlices)+2, 2), NameSuffix, '.tif'];
  imwrite(BlankImage, [OutputFolder, filesep, NewName]);
  %Copy the rest of the images
  
  if(experimentType1 || experimentType2 || experimentType3)
    slicesCounter = 1;        
    firstImage = (framesIndex-1) * NSlices(seriesIndex) * NChannels + 1 + (channelIndex - 1);
    lastImage = framesIndex * NSlices(seriesIndex) * NChannels;
    for imageIndex = firstImage:NChannels:lastImage
      if slicesCounter <= min(NSlices)
          NewName = [Prefix, '_', iIndex(numberOfFrames,3), '_z', iIndex(slicesCounter + 1, 2), NameSuffix, '.tif'];
          imwrite(LIFImages{seriesIndex}{imageIndex,1}, [OutputFolder, filesep, NewName]);
          slicesCounter = slicesCounter + 1;
      end
    end
  end
end
