% NL: Updated this drawing heavily from processLIFExportMode
% Last updated: 2020-09-03

function FrameInfo = processLSMData(Folder, D, FrameInfo, ExperimentType, ...
    Channels, ProjectionType,Prefix, OutputFolder,nuclearGUI, zslicesPadding)
    
    disp('Exporting movie file...');
    
    cleanupObj = onCleanup(@myCleanupFun);

    liveExperiment = LiveExperiment(Prefix);

    resultsFolder = liveExperiment.resultsFolder;
        
    if ~shouldExportMovieFiles
        FrameInfo = [];
    end
    moviePrecision = 'uint16';
    hisPrecision = 'uint16';
    
    %Load the reference histogram for the fake histone channel
    load('ReferenceHist.mat', 'ReferenceHist');
    
    shouldMakeMovieMat = shouldExportMovieFiles ||...
        (shouldExportNuclearProjections &&...
        ~exist([PreProcFolder,filesep, Prefix, '_movieMat.mat'],'file'));
    
    % get basic info
    NSeries = length(D);      
    
    % This chunk makes FrameInfo
    if shouldMakeMovieMat                   
      Frame_Times = []; % Store the frame information
      AllLSMImages = cell(1, NSeries);
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

        try
            StartingTime(LSMIndex) = obtainZeissStartingTime(Folder, LSMIndex, LSMMeta2, NDigits);
        catch
            StartingTime(LSMIndex) = obtainZeissStartingTime(Folder, LSMIndex, LSMMeta2, NDigits+1);
        end
        [ValueField, Frame_Times] = obtainZeissFrameTimes(LSMMeta, NSlices, LSMIndex, NPlanes, NChannels, StartingTime, Frame_Times);
        [~, FrameInfo] = createZeissFrameInfo(LSMIndex, NFrames, NSlices, FrameInfo, LSMMeta, Frame_Times, ValueField);
      end
    end
    
    if ~skipExtraction
      
        if shouldMakeMovieMat
          %Copy the data
          waitbarFigure = waitbar(0, 'Extracting LSM images');


          %Counter for number of frames
          numberOfFrames = 1;


          ySize = size(AllLSMImages{1}{1,1}, 1);
          xSize = size(AllLSMImages{1}{1,1}, 2);
          BlankImage = uint16(zeros(ySize, xSize));
  
          hisMat = zeros(ySize, xSize, sum(NFrames), hisPrecision);
        end

        if shouldExportMovieFiles
          % We need a second pass to set the correct slices count after having
          % processed all the series so we know the max(NSlices) number
          % if zPadding was indicated in the arguments, we round up to the series
          % with more z-slices (because we'll pad with blank images the other series)
%           if (zslicesPadding)
%             topZSlice = max(NSlices);
%           else
%             % if no zPadding, we round down to the series with less z-slices
%             topZSlice = min(NSlices);
%           end
          topZSlice = min(NSlices);
          
          for frameInfoIndex = 1:size(FrameInfo, 2)
            FrameInfo(frameInfoIndex).NumberSlices = topZSlice;
          end         

%           [coatChannel, histoneChannel, fiducialChannel, inputProteinChannel, FrameInfo] =...
%           LIFExportMode_interpretChannels(ExperimentType, Channels{1}, Channels{2}, Channels{3}, FrameInfo);
          
          if nuclearGUI
            warning('NL: Nuclear GUI currently not supported for LSM export.')             
          end
          
          for seriesIndex = 1:NSeries
%               for framesIndex = 1:NFrames(seriesIndex) 

%                 processMovieFrame(numberOfFrames, Prefix, OutputFolder,...
%                     AllLSMImages, framesIndex, seriesIndex, NChannels(1), NSlices, ...
%                     Channels, ProjectionType, ReferenceHist, zslicesPadding, 0);
% 
%                 numberOfFrames = numberOfFrames + 1;
%               end

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
                            
                            imwrite(AllLSMImages{seriesIndex}{imageIndex,1},...
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
                if ~nuclearGUI
                    
                    hisMat(:, :, numberOfFrames) = generateNuclearChannel(...
                        numberOfFrames, AllLSMImages,...
                        framesIndex, seriesIndex, NSlices, NChannels,ProjectionType,...
                        Channels, ReferenceHist, PreProcFolder, Prefix);
                    
                    saveNuclearProjection(hisMat, [PreProcFolder, filesep, Prefix, '-His.tif']);
                    
                end
                
                numberOfFrames = numberOfFrames + 1;
            end
          end
          
          try close(waitbarFigure); catch; end
          
          [FFPaths, FFToUse, LSMFF] = findFlatFieldInformation(Folder);
          processFlatFieldInformation(Prefix, OutputFolder, FFPaths, FFToUse, LSMFF);
        end
    end
end

