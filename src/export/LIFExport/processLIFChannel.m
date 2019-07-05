function processLIFChannel(ExperimentType, channelIndex, numberOfFrames, Prefix, BlankImage, OutputFolder,...
  LIFImages, framesIndex, seriesIndex, NChannels, NSlices, coatChannel, inputProteinChannel, zslicesPadding)
  
  experimentType1 = (strcmpi(ExperimentType,'1spot') || strcmp(ExperimentType,'2spot') || strcmp(ExperimentType,'2spot1color')) && channelIndex == coatChannel;
  experimentType2 = strcmpi(ExperimentType,'2spot2color') || strcmpi(ExperimentType, 'inputoutput');
  experimentType3 = strcmpi(ExperimentType, 'input') && sum(channelIndex == inputProteinChannel);

  % if zPadding was indicated in the arguments, we round up to the series
  % with more z-slices (because we'll pad with blank images the other series)
  if (zslicesPadding)
    topZSlice = max(NSlices);
  else
    % if no zPadding, we round down to the series with less z-slices
    topZSlice = min(NSlices);
  end
  
  % Save the blank image at the beginning of the stack
  NameSuffix = ['_ch',iIndex(channelIndex,2)];
  NewName = [Prefix, '_', iIndex(numberOfFrames,3), '_z', iIndex(1,2), NameSuffix, '.tif'];
  imwrite(BlankImage, [OutputFolder, filesep, NewName]);
  
  %Copy the rest of the images
  if(experimentType1 || experimentType2 || experimentType3)
    slicesCounter = 1;        
    firstImage = (framesIndex-1) * NSlices(seriesIndex) * NChannels + 1 + (channelIndex - 1);
    lastImage = framesIndex * NSlices(seriesIndex) * NChannels;
    for imageIndex = firstImage:NChannels:lastImage
      if slicesCounter <= topZSlice
          % if zPadding, it will process all images (because topZSlice would be max(NSlices)
          % if no zPadding, it will process images rounding down to the series with least
          % zSlices, because topZSlice would be min(NSlices)
          NewName = [Prefix, '_', iIndex(numberOfFrames,3), '_z', iIndex(slicesCounter + 1, 2), NameSuffix, '.tif'];
          imwrite(LIFImages{seriesIndex}{imageIndex,1}, [OutputFolder, filesep, NewName]);
          slicesCounter = slicesCounter + 1;
      end
    end
  end
  
  % Save as many blank images at the end of the stack are needed
  % (depending on zPadding being active or not)
  for zPaddingIndex = slicesCounter+1:topZSlice+2
    NewName = [Prefix, '_', iIndex(numberOfFrames,3), '_z', iIndex(zPaddingIndex, 2), NameSuffix, '.tif'];
    imwrite(BlankImage, [OutputFolder, filesep, NewName]);
  end
  
end
