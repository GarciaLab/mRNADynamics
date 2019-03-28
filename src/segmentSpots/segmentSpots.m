% segmentSpots(Prefix, Threshold, [Options])
%
% DESCRIPTION
% Identify and segment individual transcription spots.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% Threshold: Threshold to be used (>=). Should be kept at ~90-200 for lattice
%           light-sheet data, and at ~5-10 for confocal data (Leica SP8).
%           If left empty, then the code just generates the DoG files.
% [Options]: See below.
%
% OPTIONS
% 'displayFigures':   If you want to display plots and images.
% 'Weka': For Weka machine learning. 
%
% 'InitialFrame', N: Run the code from frame N to last frame. Defaults to first
%                frame.
%
% 'LastFrame', M:     Run the code from initial frame to frame M. Defaults to all
%                frames. It's suggested to run 5-20 frames for debugging.
%
% 'Shadows':    	 This option should be followed by 0, 1 or 2. This
%                specifies the number of requisite z-planes above and/or below the
%                brightest plane for a spot to have to pass quality control.
% 'keepPool': Don't shut down the parallel pool when the script is done
% running.
% 'nWorkers': Specify the number of workers to use during parallel
% processing
% 'noIntegralZ':  Don't establish center slice at position that maximizes raw fluo integral 
%                 across sliding 3 z-slice window.
% 'intScale': Scale up the radius of integration
% 'autoThresh': Pops up a UI to help decide on a threshhold
% 'keepProcessedData': Keeps the ProcessedData folder for the given prefix after running segment spots
% 'fit3D': Fit 3D Gaussians to all segmented spots. 
% 'skipChannel': Skips segmentation of channels inputted array (e.g. [1]
%                skips channel 1, [1, 2] skips channels 1 and 2
%
% OUTPUT
% 'Spots':  A structure array with a list of detected transcriptional loci
% in each frame and their properties.
% 'log.mat': A cell array containing logging data from the segmentation
% process. There's one row per run of segmentSpots(ML) on that particular
% dataset. 
%
% Author (contact): Armando Reimer (areimer@berkeley.edu)
% Created: 01/01/2016
% Last Updated: 8/23/2018
%
% Documented by: Armando Reimer (areimer@berkeley.edu)

function log = segmentSpots(Prefix, Threshold, varargin)

  warning('off', 'MATLAB:MKDIR:DirectoryExists');

  disp('Segmenting spots...')
  
  [displayFigures, numFrames, numShadows, intScale, nWorkers, keepPool, ...
    pool, autoThresh, initialFrame, useIntegralCenter, Weka, keepProcessedData, fit3D, skipChannel] = determineSegmentSpotsOptions(varargin);
      
  argumentErrorMessage = 'Please use filterMovie(Prefix, options) instead of segmentSpots with the argument "[]" to generate DoG images';
  try 
    if autoThresh
        Threshold = -1;
    elseif isempty(Threshold)
      error(argumentErrorMessage);
    end 

  catch 
    error(argumentErrorMessage);
  end 

 
  if nWorkers ~= 0 && ~displayFigures
    
      maxWorkers = nWorkers;

    try 
      parpool(maxWorkers); 
    catch 

      try 
        parpool; % in case there aren't enough cores on the computer
      catch 
        % parpool throws an error if there's a pool already running.
      end 

    end 

  end 

  [~, ~, ~, ~, ~, ~, ~, ExperimentType, Channel1, Channel2, ~] = readMovieDatabase(Prefix);

  [~, ProcPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

  load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat'], 'FrameInfo');

  ProcessedDataFolder = [ProcPath, filesep, Prefix, '_'];
  DogOutputFolder = [ProcessedDataFolder, filesep, 'dogs'];
  mkdir(DogOutputFolder)

  microscope = FrameInfo(1).FileMode;

  zSize = FrameInfo(1).NumberSlices + 2;

  if numFrames == 0
    numFrames = length(FrameInfo);
  end 

  nCh = 1;

  if strcmpi(ExperimentType, '2spot2color')
    nCh = 2;
  end 

  [ffim, doFF] = loadSegmentSpotsFlatField(PreProcPath, Prefix, FrameInfo);
  clear rawdir;

  % The spot finding algorithm first segments the image into regions that are
  % above the threshold. Then, it finds global maxima within these regions by searching in a region "neighborhood"
  % within the regions.

  pixelSize = FrameInfo(1).PixelSize * 1000; %nm
  neighborhood = round(1300 / pixelSize); %nm
  snippet_size = 2 * (floor(1300 / (2 * pixelSize))) + 1; % nm. note that this is forced to be odd

  all_frames = cell(numFrames, zSize);
  close all;

  coatChannel = getCoatChannel(ExperimentType, Channel1, Channel2);
  Spots = [];
  falsePositives = 0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Segment transcriptional loci
      
  for channelIndex = 1:nCh
      if ismember(channelIndex, skipChannel)
          continue
      end
      
    tic;
      
    all_frames = segmentTranscriptionalLoci(ExperimentType, coatChannel, channelIndex, all_frames, initialFrame, numFrames, zSize, ...
      PreProcPath, Prefix, DogOutputFolder, displayFigures, pool, doFF, ffim, Threshold(channelIndex), neighborhood, ...
      snippet_size, pixelSize, microscope, intScale, Weka);

    close all;

    % Create a useful structure that can be fed into pipeline
    [Particles, fields] = saveParticleInformation(numFrames, all_frames, zSize, useIntegralCenter);

    if ~isempty(fields)
       
        [neighborhood, Particles] = segmentSpotsZTracking(pixelSize, numFrames, Particles, fields); 

        [Particles, falsePositives] = findBrightestZ(Particles,numShadows, useIntegralCenter, 0);

        %Create a final Spots structure to be fed into TrackmRNADynamics
        Spots{channelIndex} = createSpotsStructure(Particles, numFrames, 1);
 
    end
    
    timeElapsed = toc;
    disp(['Elapsed time: ', num2str(timeElapsed / 60), ' min'])
    try %#ok<TRYNC>
        log = logSegmentSpots(DropboxFolder, Prefix, timeElapsed, [], numFrames, Spots, falsePositives, Threshold, channelIndex);
        display(log);
    end
    
  end
  
    %If we only have one channel, then convert Spots to a
    %standard structure.
  if nCh == 1 && iscell(Spots)
      Spots = Spots{1}; %#ok<NASGU>
  end 
  
  mkdir([DropboxFolder, filesep, Prefix]);
  save([DropboxFolder, filesep, Prefix, filesep, 'Spots.mat'], 'Spots', '-v7.3');


%   t = toc;
%   disp(['Elapsed time: ', num2str(t / 60), ' min'])
% 
%   log = logSegmentSpots(DropboxFolder, Prefix, t, numFrames, Spots, falsePositives, Threshold);
%   display(log);

  if ~keepProcessedData
    deleteProcessedDataFolder(ProcessedDataFolder, Prefix);
  else
    disp('keepProcessedData parameter sent. ProcessedData folder will not be removed.');
  end

  if ~keepPool

    try  %#ok<TRYNC>
      poolobj = gcp('nocreate');
      delete(poolobj);
    end 

  end 

    if fit3D
        disp('Fitting 3D Gaussians...')
        fit3DGaussiansToAllSpots(Prefix, 'segmentSpots', Spots);
        disp('3D Gaussian fitting completed.')
    end
  
end 
