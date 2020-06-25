function FrameInfo = processLSMData(Folder, D, FrameInfo, ExperimentType, ...
    Channels, ProjectionType,Prefix, OutputFolder,nuclearGUI, zslicesPadding)
  % What type of experiment do we have?

    NSeries = length(D);
    Frame_Times = []; % Store the frame information
    AllLSMImages = cell(1, NSeries);
    
    % preprocess data
    waitbarFigure = waitbar(0, 'Extracting LSM images');
    load('ReferenceHist.mat')
    for LSMIndex = 1:NSeries
      waitbar(LSMIndex / NSeries, waitbarFigure);
      
      %Load the file using BioFormats
      LSMImages = bfopen([Folder, filesep, D(LSMIndex).name]);
      % Extract the metadata for each series
      LSMMeta = LSMImages{:, 4}; % OME Metadata
      LSMMeta2 = LSMImages{:, 2}; % Original Metadata
      AllLSMImages{LSMIndex} = LSMImages{1};

      % Figure out the number of slices in each series
      NSlices(LSMIndex) = str2double(LSMMeta.getPixelsSizeZ(0));
      % Number of channels
      NChannels(LSMIndex) = LSMMeta.getChannelCount(0);
      % Total number of planes acquired
      NPlanes(LSMIndex) = LSMMeta.getPlaneCount(0);
      % Finally, use this information to determine the number of frames in
      % each series
      NFrames(LSMIndex) = NPlanes(LSMIndex) / NSlices(LSMIndex) / NChannels(LSMIndex);

      % Check that the acquisition wasn't stopped before the end of a
      % cycle. If that is the case, some images in the last frame will
      % be blank and we need to remove them.
      if sum(sum(LSMImages{1}{end,1})) == 0
        % Reduce the number of frames by one
        NFrames(LSMIndex) = NFrames(LSMIndex) - 1;
        % Reduce the number of planes by NChannels*NSlices
        NPlanes(LSMIndex) = NPlanes(LSMIndex) - NChannels(LSMIndex) * NSlices(LSMIndex);
      end

      NDigits = getNDigits(NFrames, LSMIndex);

      StartingTime(LSMIndex) = obtainZeissStartingTime(Folder, LSMIndex, LSMMeta2, NDigits);
      [ValueField, Frame_Times] = obtainZeissFrameTimes(LSMMeta, NSlices, LSMIndex, NPlanes, NChannels, StartingTime, Frame_Times);
      [~, FrameInfo] = createZeissFrameInfo(LSMIndex, NFrames, NSlices, FrameInfo, LSMMeta, Frame_Times, ValueField);
    end
    
    
    % We need a second pass to set the correct slices count after having
    % processed all the series so we know the max(NSlices) number
    % if zPadding was indicated in the arguments, we round up to the series
    % with more z-slices (because we'll pad with blank images the other series)
    if (zslicesPadding)
      topZSlice = max(NSlices);
    else
      % if no zPadding, we round down to the series with less z-slices
      topZSlice = min(NSlices);
    end
    
    for frameInfoIndex = 1:size(FrameInfo, 2)
      FrameInfo(frameInfoIndex).NumberSlices = topZSlice;
    end
  
    close(waitbarFigure);
    
    [coatChannel, histoneChannel, fiducialChannel, inputProteinChannel, FrameInfo] =...
    LIFExportMode_interpretChannels(ExperimentType, Channels{1}, Channels{2}, Channels{3}, FrameInfo);

    numberOfFrames = 1;
    %at the moment, this doesn't work. it would be nice to correctly
    %implement it in the future. -AR 6/25/2020
%     if nuclearGUI
%         [~, ProjectionType] = chooseNuclearChannels(...
%         AllLSMImages, NSeries, NSlices, NChannels(1), NFrames, ProjectionType, Channel1, Channel2, ...
%         Channel3, ReferenceHist);
%     end
    for seriesIndex = 1:NSeries
        for framesIndex = 1:NFrames(seriesIndex) 

          processMovieFrame(numberOfFrames, Prefix, OutputFolder,...
              AllLSMImages, framesIndex, seriesIndex, NChannels(1), NSlices, ...
              Channels, ProjectionType, ReferenceHist, zslicesPadding, 0);
          
          numberOfFrames = numberOfFrames + 1;
        end
    end

    [FFPaths, FFToUse, LSMFF] = findFlatFieldInformation(Folder);
    processFlatFieldInformation(Prefix, OutputFolder, FFPaths, FFToUse, LSMFF);

end

