% segmentSpots(Prefix, Threshold, [Options])
%
% DESCRIPTION
% Identify and segment individual transcription spots.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% Threshold: Threshold to be used. Should be kept at ~90-200 for lattice
%           light-sheet data, and at ~5-10 for confocal data (Leica SP8).
%           If left empty, then the code just generates the DoG files.
% [Options]: See below.
%
% OPTIONS
% 'displayFigures':   If you want to display plots and images.
%
% 'Frames', N:Run the code from frame 1 to frame N. Defaults to all
%                frames. It's suggested to run 5-20 frames for debugging.
%
% 'Shadows':    	 This option should be followed by 0, 1 or 2. This
%                specifies the number of requisite z-planes above and/or below the
%                brightest plane for a spot to have to pass quality control.
% 'noPool':     Does not start and use a parallel pool.
% 'keepPool': Don't shut down the parallel pool when the script is done
% running.
% 'highPrecision': Uses higher precision filtering for segmentation
% 'intScale': Scale up the radius of integration
% 'customFilters': Choose which filter to use to segment the image. Name
%                  should be a string, followed by a cell with your filter
%                  or filters
%                 ex. segmentSpots(Prefix,[],'customFilter', 'Structure_largest', {1, 8})
%           Filter Options:
%               'Gaussian_blur''Median'
%               'Edges''Maximum'
%               'Laplacian''Minimum'
%               'Mean''Std'
%               'Hessian_largest''Hessian_smallest'
%               [DEFAULT] 'Difference_of_Gaussian' (2 sigmas) [DEFAULT]
%               'Structure_largest' (2 sigmas)
%               'Structure_smallest' (2 sigmas)
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
% Last Updated: 8/17/2018
%
% Documented by: Armando Reimer (areimer@berkeley.edu)

function log = segmentSpots(Prefix, Threshold, varargin)

  warning('off', 'MATLAB:MKDIR:DirectoryExists');

  [displayFigures, trackSpots, numFrames, numShadows, customFilter, highPrecision, filterType, ...
  intScale, nWorkers, keepPool, pool] = determineSegmentSpotsOptions(varargin);

  % If no threshold was specified, then just generate the DoG images
  justDoG = 0;

  try 

    if isempty(Threshold)
      justDoG = 1;
    end 

  catch 
    error('Please pass the argument "[]" to generate DoG images')
  end 

  % Start timer
  tic;

  if pool
    maxWorkers = nWorkers;

    try 
      parpool(maxWorkers); % 6 is the number of cores the Garcia lab server can reasonably handle per user at present.
    catch 

      try 
        parpool; % in case there aren't enough cores on the computer
      catch 
        % parpool throws an error if there's a pool already running.
      end 

    end 

  end 

  [~, ~, ~, ~, ~, ~, ~, ExperimentType, Channel1, Channel2, ~] = readMovieDatabase(Prefix);

  [~, FISHPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

  load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat']);

  % TODO harrypotel Should be renamed to DoGsFolder or something clearer.
  DogOutputFolder = [FISHPath, filesep, Prefix, '_', filesep, 'dogs'];
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
  close all force;

  coatChannel = determineSegmentSpotsCoatChannel(ExperimentType, Channel1, Channel2);
  Spots = [];
  sigmas = [];
  falsePositives = 0;

  if justDoG
    % Generate difference of Gaussian images if no threshold was given.
    [sigmas] = generateDifferenceOfGaussianImages(DogOutputFolder, pixelSize, customFilter, nCh, ExperimentType, ...
      coatChannel, numFrames, displayFigures, zSize, PreProcPath, Prefix, filterType, highPrecision);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Segment transcriptional loci
  else 
      
    for channelIndex = 1:nCh
      all_frames = segmentTranscriptionalLoci(ExperimentType, coatChannel, channelIndex, all_frames, numFrames, zSize, ...
        PreProcPath, Prefix, DogOutputFolder, displayFigures, pool, doFF, ffim, Threshold, neighborhood, ...
        snippet_size, pixelSize, microscope, intScale);

      close all force;

      % Create a useful structure that can be fed into pipeline
      [Particles, fields] = saveParticleInformation(numFrames, all_frames, zSize);

      [neighborhood, Particles] = segmentSpotsZTracking(pixelSize, numFrames, Particles, fields); %#ok<ASGLU>

      [Particles, falsePositives] = findBrightestZ(numShadows, 0, 0);

      %Create a final Spots structure to be fed into TrackmRNADynamics
      Spots = createSpotsStructure(Particles, numFrames, channelIndex);

      %If we only have one channel, then convert Spots to a
      %standard structure.

      if nCh == 1
        Spots = Spots{1};
      end 

      mkdir([DropboxFolder, filesep, Prefix]);
      save([DropboxFolder, filesep, Prefix, filesep, 'Spots.mat'], 'Spots', '-v7.3');
    end 

  end 

  t = toc;
  disp(['Elapsed time: ', num2str(t / 60), ' min'])

  log = logSegmentSpots(DropboxFolder, Prefix, t, numFrames, justDoG, Spots, filterType, sigmas, falsePositives, Threshold);

  if ~ keepPool

    try 
      poolobj = gcp('nocreate');
      delete(poolobj);
    catch 
    end 

  end 

end 
