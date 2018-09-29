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
% 'keepPool': Don't shut down the parallel pool when the script is done
% running.
% 'nWorkers': Specify the number of workers to use during parallel
% processing
% 'IntegralZ':  Establish center slice at position that maximizes raw fluo integral 
%               across sliding 3 z-slice window.
% 'intScale': Scale up the radius of integration
% 'autoThresh': Pops up a UI to help decide on a threshhold
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

  [displayFigures, numFrames, numShadows, intScale, nWorkers, keepPool, ...
      pool, autoThresh] = determineSegmentSpotsOptions(varargin);

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

  [~, FISHPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

  load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat']);

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
  use_integral_center = 1; %for z-tracking

  all_frames = cell(numFrames, zSize);
  close all;

  coatChannel = getCoatChannel(ExperimentType, Channel1, Channel2);
  Spots = [];
  falsePositives = 0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Segment transcriptional loci
      
  for channelIndex = 1:nCh
      
    tic;
      
    all_frames = segmentTranscriptionalLoci(ExperimentType, coatChannel, channelIndex, all_frames, numFrames, zSize, ...
      PreProcPath, Prefix, DogOutputFolder, displayFigures, pool, doFF, ffim, Threshold, neighborhood, ...
      snippet_size, pixelSize, microscope, intScale);

    close all;

    % Create a useful structure that can be fed into pipeline
    [Particles, fields] = saveParticleInformation(numFrames, all_frames, zSize, use_integral_center);

    if ~isempty(fields)
       
        [neighborhood, Particles] = segmentSpotsZTracking(pixelSize, numFrames, Particles, fields); 

        [Particles, falsePositives] = findBrightestZ(Particles,numShadows, 0, 0);

        %Create a final Spots structure to be fed into TrackmRNADynamics
        Spots{channelIndex} = createSpotsStructure(Particles, numFrames);
 
    end
    
    t = toc;
    disp(['Elapsed time: ', num2str(t / 60), ' min'])
    try %#ok<TRYNC>
        log = logSegmentSpots(DropboxFolder, Prefix, t, numFrames, Spots, falsePositives, Threshold, channelIndex);
        display(log);
    end
    
  end
  
    %If we only have one channel, then convert Spots to a
    %standard structure.
  if nCh == 1
      Spots = Spots{1}; %#ok<NASGU>
  end 
  
  mkdir([DropboxFolder, filesep, Prefix]);
  save([DropboxFolder, filesep, Prefix, filesep, 'Spots.mat'], 'Spots', '-v7.3');


%   t = toc;
%   disp(['Elapsed time: ', num2str(t / 60), ' min'])
% 
%   log = logSegmentSpots(DropboxFolder, Prefix, t, numFrames, Spots, falsePositives, Threshold);
%   display(log);

  if ~keepPool

    try  %#ok<TRYNC>
      poolobj = gcp('nocreate');
      delete(poolobj);
    end 

  end 

end 
