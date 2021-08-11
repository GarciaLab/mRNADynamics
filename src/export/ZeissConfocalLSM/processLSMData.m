% NL: Updated this drawing heavily from processLIFExportMode
% Last updated: 2020-09-03

function FrameInfo = processLSMData(Folder, RawDataFiles, FrameInfo,...
    Channels, ProjectionType, Prefix, OutputFolder,nuclearGUI,...
    skipExtraction,skipNuclearProjection,zslicesPadding)
    
    disp('Exporting movie file...');
    
    cleanupObj = onCleanup(@myCleanupFun);
    
    %Load the reference histogram for the fake histone channel
    load('ReferenceHist.mat', 'ReferenceHist');    
    
    % initialize FrameInfo
    FrameInfo = [];
        
    % get basic info
    NSeries = length(RawDataFiles);      
        
    if ~skipExtraction
      % This chunk makes FrameInfo                     
      [FrameInfo,AllLSMImages,NSlices, ~, NFrames,~,NChannels] ...
        = getZeissFrameInfo(RawDataFiles,NSeries,FrameInfo,zslicesPadding);
      
      % pop channels up a level
      if iscell(Channels{1})
          Channels = [Channels{:}];
      end
      NChannelsAlt = find(~cellfun(@isempty,Channels),1,'last');
      if NChannelsAlt~=NChannels
          warning(['issue with NChannels metadata extraction. Reseting to ' num2str(NChannelsAlt)])
          NChannels = NChannelsAlt;
      end
      
      % save FrameInfo
      liveExperiment = LiveExperiment(Prefix);
      save([liveExperiment.resultsFolder,filesep,'FrameInfo.mat'], 'FrameInfo');
      
      % this function exports tif z stacks
      exportTifStacks(AllLSMImages, 'LSM', NChannels, NFrames, NSlices, Prefix, ...
          ...moviePrecision, hisPrecision, 
          nuclearGUI, ProjectionType, Channels, ReferenceHist, ...
          skipNuclearProjection,zslicesPadding);
        
      % Look for flat field images
      [FFPaths, FFToUse, LSMFF] = findFlatFieldInformation(Folder);
      
      % Proceed accordingly
      processFlatFieldInformation(Prefix, OutputFolder, FFPaths, FFToUse, LSMFF);

      imarisFolder = [RawDataFiles(1).folder, filesep, 'ImarisResult'];
      if isfolder(imarisFolder)
          disp('Imaris results folder detected. Will get nuclear tracking info from Imaris.')
          
          DropboxFolder = liveExperiment.userResultsFolder;
          ellipsesFile = [DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'];
          schnitzcellsFile = [DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat']; 

          imarisStatisticsFolder = [imarisFolder, filesep, liveExperiment.experimentName, '_Statistics'];
          positionsFile = [imarisStatisticsFolder, filesep, liveExperiment.experimentName, '_Position.csv'];

          % multi-file version of imaris data where there's a statistics subfolder with a positions file
          % for each movie file instead of a single csv file for the whole movie
          multiFileStats = dir([imarisFolder, filesep, '**/*Out_Position.csv']);
          
          if isfile(positionsFile)
            disp('Parsing Imaris position file...')

            

            % Liz says: The entire image stack is 200 pixels in x, 444 in y and 100 in z (again will change in future datasets).
            % The x and y resolutions is 0.234 microns per pixel and z is 0.5 microns per pixel, but this will change in future datasets.
            % I looked up this information in the Zeiss software, but hopefully it's available in the metadata. 
            % JP: it is, in FrameInfo:

            % LinesPerFrame: 444
            % PixelsPerLine: 200
            % NumberSlices: 100
            % FileMode: 'LSMExport'
            % PixelSize: 0.2341
            % ZStep: 0.5000
            % Time: 1.4797e+04

            [schnitzcells, Ellipses] = readimariscsv(0, imarisStatisticsFolder, positionsFile, FrameInfo(1).PixelSize, FrameInfo(1).PixelSize, FrameInfo(1).ZStep);

            save2(ellipsesFile, Ellipses); 
            save2(schnitzcellsFile, schnitzcells); 

            disp('Saved imaris-based Ellipses and schnitzcells files in DynamicsResults folder.');
          elseif ~isempty(multiFileStats)
            disp('Parsing mult-file Imaris statistic data...')

            % because each Imaris file corresponds to a different movie sub-file,
            % we need to offset frame numbers to match overall frame numbering
            % for example file 1 contains frames 1 and 2, and file 2 contains also frames 1 and 2,
            % but for us it's file 1 = frame 1 and 2, file 2 = frames 3 and 4, and we schnitzcells and Ellipses
            % to reflect that numbering
            frameNumberOffset = 0;
            
            for positionFileIndex = 1:size(multiFileStats)
              positionSubFile = multiFileStats(positionFileIndex);
              fprintf('Parsing position file at %s\n', positionSubFile.name);
              
              % each csv processing will tell us the last frame number it referenced
              % so we can use that to offset the frame numbers in next file
              [partialSchnitzcells, partialEllipses, lastFrame] = readimariscsv(frameNumberOffset, positionSubFile.folder, [positionSubFile.folder, filesep, positionSubFile.name], FrameInfo(1).PixelSize, FrameInfo(1).PixelSize, FrameInfo(1).ZStep);

              % append partial info from schnitcells and ellipses from sub file
              % to the general structs
              if (positionFileIndex == 1)
                % it's the first file, create general variables
                Ellipses = partialEllipses;
                schnitzcells = partialSchnitzcells;
              else
                % it's not the first file, variables already exist, so we concatenate to them
                schnitzcells = [schnitzcells; partialSchnitzcells];
                Ellipses = [Ellipses; partialEllipses];
              end

              frameNumberOffset = lastFrame;              
              fprintf('Last frame number was %d\n', frameNumberOffset);

            end

            save2(ellipsesFile, Ellipses); 
            save2(schnitzcellsFile, schnitzcells); 

          else
            error(['Statistics files not found in Imaris folder.', positionsFile]);
          end
          

      end

      if nuclearGUI

        chooseAnaphaseFrames(...
            Prefix, 'ProjectionType', ProjectionType,...
            'ReferenceHist', ReferenceHist);

      end
    end
end

