function exportTifStacks(FrameInfo,AllImages, imagingModality, NChannels, ...
  Prefix, moviePrecision, hisPrecision, nuclearGUI, ProjectionType, Channels, ReferenceHist)

  % extract metadata from FrameInfo
  NSlices = FrameInfo(1).NumberSlices;
  NSeries = length(AllImages);
  NFrames = length(FrameInfo);
  
  liveExperiment = LiveExperiment(Prefix);
  PreProcFolder = liveExperiment.preFolder;
% if shouldMakeMovieMat
  %Copy the data
  waitbarFigure = waitbar(0, ['Extracting ' imagingModality ' images']);


  %Counter for number of frames
  numberOfFrames = 1;


  ySize = size(AllImages{1}{1,1}, 1);
  xSize = size(AllImages{1}{1,1}, 2);
  BlankImage = zeros(ySize, xSize, moviePrecision);

  hisMat = zeros(ySize, xSize, sum(NFrames), hisPrecision);
% end


% if shouldExportMovieFiles
        
  topZSlice = min(NSlices);

  for seriesIndex = 1:NSeries
      waitbar(seriesIndex/NSeries, waitbarFigure)
      for framesIndex = 1:NFrames(seriesIndex)

          for channelIndex = 1:NChannels

              NameSuffix = ['_ch',iIndex(channelIndex,2)];

              NewName = [Prefix, '_', iIndex(numberOfFrames,3),...
                  NameSuffix, '.tif'];

              imwrite(BlankImage, [PreProcFolder, filesep, NewName]);
              %
              %Copy the rest of the images
              slicesCounter = 1;
              firstImageIndex = (framesIndex-1) * NSlices(seriesIndex) * NChannels +...
                  1 + (channelIndex - 1);
              lastImageIndex = framesIndex * NSlices(seriesIndex) * NChannels;
              if firstImageIndex == lastImageIndex
                  firstImageIndex = 1;
                  lastImageIndex = 1;
              end
              for imageIndex = firstImageIndex:NChannels:lastImageIndex
                  if slicesCounter <= topZSlice
                      % if zPadding, it will process all images (because topZSlice would be max(NSlices)
                      % if no zPadding, it will process images rounding down to the series with least
                      % zSlices, because topZSlice would be min(NSlices)

                      imwrite(AllImages{seriesIndex}{imageIndex,1},...
                          [PreProcFolder, filesep, NewName], 'WriteMode', 'append');
                      slicesCounter = slicesCounter + 1;

                  end
              end

              %Save as many blank images at the end of the stack are needed
              %(depending on zPadding being active or not)
              for zPaddingIndex = slicesCounter+1:topZSlice+2
                  imwrite(BlankImage, [PreProcFolder, filesep, NewName], 'WriteMode', 'append');
              end
          end



          %Now create nuclear projection movies
          if ~nuclearGUI || strcmpi(imagingModality,'LSM')

              hisMat(:, :, numberOfFrames) = generateNuclearChannel(...
                  numberOfFrames, AllImages,...
                  framesIndex, seriesIndex, NSlices, NChannels,ProjectionType,...
                  Channels, ReferenceHist, PreProcFolder, Prefix);

              saveNuclearProjection(hisMat, [PreProcFolder, filesep, Prefix, '-His.tif']);

          end

          numberOfFrames = numberOfFrames + 1;
      end
  end


% end