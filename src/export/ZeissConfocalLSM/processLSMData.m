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
          disp('Imaris results folder detected. Will get spots info from Imaris.')
          
          positionsFile = [imarisFolder, filesep, liveExperiment.experimentName, '_Statistics', filesep, liveExperiment.experimentName, '_Position_clean.csv'];
          if isfile(positionsFile)
            disp('Parsing Imaris position file...')
            
            imarisPositions = csv2cell(positionsFile, 'fromfile');
            
          else
            error(['Imaris folder does not contain positions file at ', positionsFile]);
          end
          
      end
      % TO-DO: extract metadata for Imaris processing. Liz says:
      %  The x and y resolutions is 0.234 microns per pixel and z is 0.5 microns per pixel, but this will change in future datasets.
      % The entire image stack is 200 pixels in x, 444 in y and 100 in z (again will change in future datasets).
      % I looked up this information in the Zeiss software, but hopefully it's available in the metadata. 
      % JP: it is, in FrameInfo:

      % LinesPerFrame: 444
      % PixelsPerLine: 200
      % NumberSlices: 100
      % FileMode: 'LSMExport'
      % PixelSize: 0.2341
      % ZStep: 0.5000
      % Time: 1.4797e+04


      % TO-DO: check if folder ./imarisResults exists in RawDynamicsData folder, and if so, look for
      % subfolder named [experimentName]_Statistics and file [experimentName]_Position.csv
      % Inside that file, columns are as follows:
      % Position X  Position Y  Position Z  Unit  Category  Collection  Time  TrackID ID
      % Header row is row #4 (the first 3 rows are title and blank)
      % To-Do: check if we can clean it up prior to reading into a cell array
      % Time is the frame number. TrackID seems to be the cross-frame nuclei ID, while "ID" is just a unique ID for each row.

      % According to our "Getting Started" doc, segmentSpots comes first and TrackNuclei later. Will Imaris replace both?
      % segmentSpots output according to documentation is:
      % OUTPUT
      % 'Spots':  A structure array with a list of detected transcriptional loci in each frame and their properties.
      % (a log also, probably with dont need that)

      % Example of Spots.mat data from actual CZI data:

      % FixedAreaIntensity: [7.2018e+04 6.0188e+04 1.0406e+05 8.5536e+04]
      % xFit: [71.3392 71.5677 73.1941 73.0870]
      % yFit: [263.1056 262.7039 260.9781 261.0533]
      % Offset: [0.1000 0.1000 1.9703e+03 2.5219e+03]
      % xFitWidth: [1.9029 1.6573 2.1945 1.7687]
      % yFitWidth: [1.3755 1.3648 1.6780 1.6502]
      % yDoG: [262 262 260 260]
      % xDoG: [72 72 72 72]
      % GaussianIntensity: [3.2948e+05 3.1675e+05 1.6708e+05 165120]
      % CentralIntensity: [5876 5248 5964 6676]
      % DOGIntensity: [50604 49539 51858 52325]
      % ConfidenceIntervals: [9×8 single]
      % gaussParams: {[1×9 single]  [1×9 single]  [1×9 single]  [1×9 single]}
      % dogFixedAreaIntensity: [2459599 2436644 3165465 3307221]
      % intArea: 109
      % z: [4 5 79 80]
      % FixedAreaIntensity3: 1.3221e+05
      % brightestZ: 4
      % snippet_size: 5
      % Approved: [0 0 0 0]

      % JP: I'm guessing we'll need to process all of Imaris CSV to replace this data. Or. we could get whatever data Imaris outpus, and try to use it to make segmentSpots run easier.
      
      % TO-DO: compare what TrackNuclei outputs to what we have in this CSV and code the "adapter" from Imaris to "TrackNuclei output"
      % according to TarckNuclei documentation:
      % OUTPUT.
      % '*_lin.mat' : Nuclei with lineages
      % 'Ellipses.mat' : Just nuclei

      

      % Process Imaris and the check manually with CheckNucleiSegmentation?
      
      if nuclearGUI

        chooseAnaphaseFrames(...
            Prefix, 'ProjectionType', ProjectionType,...
            'ReferenceHist', ReferenceHist);

      end
    end
end

