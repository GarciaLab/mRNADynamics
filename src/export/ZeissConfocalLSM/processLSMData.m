% NL: Updated this drawing heavily from processLIFExportMode
% Last updated: 2020-09-03

function FrameInfo = processLSMData(Folder, RawDataFiles, FrameInfo,...
    Channels, ProjectionType, Prefix, OutputFolder,nuclearGUI,...
    skipExtraction,skipNuclearProjection,zslicesPadding)
    
    disp('Exporting movie file...');
    
    cleanupObj = onCleanup(@myCleanupFun);

    moviePrecision = 'uint16';
    hisPrecision = 'uint16';
    
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
      
      % save FrameInfo
      liveExperiment = LiveExperiment(Prefix);
      save([liveExperiment.resultsFolder,filesep,'FrameInfo.mat'], 'FrameInfo')
      
      % this function exports tif z stacks
      exportTifStacks(AllLSMImages, 'LSM', NChannels, NFrames, NSlices, Prefix, ...
          moviePrecision, hisPrecision, nuclearGUI, ProjectionType, Channels, ReferenceHist, ...
          skipNuclearProjection,zslicesPadding)           
      
      % Look for flat field images
      [FFPaths, FFToUse, LSMFF] = findFlatFieldInformation(Folder);
      % Proceed accordingly
      processFlatFieldInformation(Prefix, OutputFolder, FFPaths, FFToUse, LSMFF);

      imarisFolder = [RawDataFiles.folder, filesep, 'ImarisResult'];
      if isfolder(imarisFolder)
          disp('Imaris results folder detected. Will get nuclear tracking info from Imaris.')
          
          imarisStatisticsFolder = [imarisFolder, filesep, liveExperiment.experimentName, '_Statistics'];
          positionsFile = [imarisStatisticsFolder, filesep, liveExperiment.experimentName, '_Position.csv'];
          if isfile(positionsFile)
            disp('Parsing Imaris position file...')

            DropboxFolder = liveExperiment.userResultsFolder;
            ellipsesFile = [DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'];
            schnitzcellsFile = [DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat']; 

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

            [schnitzcells, Ellipses] = readimariscsv(imarisStatisticsFolder, positionsFile, FrameInfo(1).PixelSize, FrameInfo(1).PixelSize, FrameInfo(1).ZStep);


            save2(ellipsesFile, Ellipses); 
            save2(schnitzcellsFile, schnitzcells); 

            disp('Saved imaris-based Ellipses and schnitzcells files in DynamicsResults folder.');
          else
            error(['Imaris folder does not contain positions file at ', positionsFile]);
          end
          
      end

      if nuclearGUI

        chooseAnaphaseFrames(...
            Prefix, 'ProjectionType', ProjectionType,...
            'ReferenceHist', ReferenceHist);

      end
    end
end

