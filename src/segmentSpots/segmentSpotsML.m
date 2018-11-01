% segmentSpotsML(Prefix, Threshold, [Options])
%
% DESCRIPTION
% Identify and segment individual transcription Spots.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% Threshold: Threshold to be used. Should be kept at 5000 always.
%
% [Options]: See below.
%
% OPTIONS
% 'displayFigures':   If you want to display plots and images.
%
% 'InitialFrame', N: Run the code from frame N to last frame. Defaults to first
%                frame.
%
% 'LastFrame', M:     Run the code from initial frame to frame M. Defaults to all
%                frames. It's suggested to run 5-20 frames for debugging.
%
% 'Shadows':    	 This option should be followed by 0, 1 or 2. This
%                specifies the number of requisite z-planes above or below the
%                brightest plane for a spot to have to pass quality control.
%
% 'IntegralZ':  Establish center slice at position that maximizes raw fluo integral
%               across sliding 3 z-slice window.
%
% 'intScale': Scale up the radius of integration
%
% 'keepPool': Don't shut down the parallel pool when the script is done
% running.
%
% 'nWorkers': Specify the number of workers to use during parallel
% processing
%
% OUTPUT
% 'Spots':  A structure array with a list of detected transcriptional loci
% in each frame and their properties.
%
% Author (contact): Armando Reimer (areimer@berkeley.edu)
% Created: 01/01/2016
% Last Updated: 12/31/2017
%
% Documented by: Armando Reimer (areimer@berkeley.edu)
function segmentSpotsML(Prefix, Threshold, varargin)

  if isempty(Threshold)
    argumentErrorMessage = ['Please use filterMovie(Prefix, ''Weka'', options) instead ', ...
                                'of segmentSpots with the argument "[]" to generate DoG images'];
    error(argumentErrorMessage);
  end

  [displayFigures, numFrames, numShadows, intScale, nWorkers, keepPool, ~, ~, ...
     initialFrame, ~] = determineSegmentSpotsOptions(varargin);
  useIntegralCenter = 0;

  %% Start timer
  tic;

  [~, ~, ~, ~, ~, ~, ~, ExperimentType, Channel1, Channel2, ~] = readMovieDatabase(Prefix);
  [~, FISHPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

  load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat']);

  microscope = FrameInfo(1).FileMode;
  zSize = FrameInfo(1).NumberSlices + 2;

  if numFrames == 0
    numFrames = length(FrameInfo);
  end

  dogsFolder = [FISHPath, filesep, Prefix, '_', filesep, 'dogs'];
  mkdir(dogsFolder)

  nCh = 1;

  if strcmpi(ExperimentType, '2spot2color')
    nCh = 2;
  end

  coatChannel = getCoatChannel(ExperimentType, Channel1, Channel2);

  %Load and apply flat-field correction
  doFF = 1;

  try
    ffim = imread([PreProcPath, filesep, Prefix, filesep, Prefix, '_FF.tif']);
    ffim = CPsmooth(ffim, 'Gaussian Filter', 256, 0); %large feature scale-space representation
    ffim = double(ffim / max(max(ffim)));
  catch
    ffim = ones(FrameInfo(1).LinesPerFrame,FrameInfo(1).PixelsPerLine);
    warning('Will not apply flat field correction');
    doFF = 0;
  end

  clear rawdir;
  % The spot finding algorithm first segments the image into regions that are
  % above the threshold. Then, it finds global maxima within these regions by searching in a region "neighborhood"
  % within the regions.

  pixelSize = FrameInfo(1).PixelSize * 1000; %nm
  neighborhood = round(1300 / pixelSize); %nm
  neighborhoodZ = neighborhood; %nm
  snippet_size = 2 * (floor(1300 / (2 * pixelSize))) + 1; % nm. note that this is forced to be odd

  allFrames = cell(numFrames, zSize);
  close all force;

  if nWorkers > 0
    maxWorkers = nWorkers;
    p = gcp('nocreate');

    if isempty(p)

      try
        parpool(maxWorkers);
      catch
        parpool;
      end

    elseif p.NumWorkers > maxWorkers
      delete(gcp('nocreate')); % if pool with too many workers, delete and restart

      try
        parpool(maxWorkers);
      catch
        parpool;
      end

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Segment transcriptional loci

  for channelIndex = 1:nCh

    if strcmpi(ExperimentType, 'inputoutput')
      nameSuffix = ['_ch', iIndex(coatChannel, 2)];
    else
      nameSuffix = ['_ch', iIndex(channelIndex, 2)];
    end

    segmentSpotsWaitbar = waitbar(0, 'Segmenting Spots');

    allFrames = processFramesSegmentSpotsML(Prefix, Threshold, PreProcPath, dogsFolder, channelIndex, nameSuffix, ...
      initialFrame, numFrames, zSize, doFF, ffim, displayFigures, nWorkers, segmentSpotsWaitbar, allFrames, ...
      snippet_size, neighborhood, microscope, intScale, pixelSize);

    close(segmentSpotsWaitbar);
    close all force;

    %Create a useful structure that can be fed into pipeline
    [Particles, particleFields] = createParticlesStructure(initialFrame, numFrames, allFrames, zSize, snippet_size, useIntegralCenter);

    %z-tracking
    Particles = findZColumns(Particles, particleFields, initialFrame, numFrames, neighborhoodZ);

    %pick the brightest z-slice
    [Particles, falsePositives] = findBrightestZ(Particles, numShadows, useIntegralCenter, 0);

    %Create a final Spots structure to be fed into TrackmRNADynamics
    Spots{channelIndex} = createSpotsStructure(Particles, numFrames, initialFrame);

    %If we only have one channel, then convert Spots{channelIndex} to a
    %standard structure.
    if nCh == 1
      Spots = Spots{1};
    end

    mkdir([DropboxFolder, filesep, Prefix]);
    save([DropboxFolder, filesep, Prefix, filesep, 'Spots.mat'], 'Spots', '-v7.3');

    t = toc;
    disp(['Elapsed time: ', num2str(t / 60), ' min'])

    try %#ok<TRYNC>
      log = logSegmentSpots(DropboxFolder, Prefix, t, initialFrame, numFrames, Spots, falsePositives, Threshold, channelIndex);
      display(log);
    end
  
  end

  if ~ keepPool

    try
      poolobj = gcp('nocreate');
      delete(poolobj);
    catch
      %fails if the parallel pool has timed out.
    end

  end

end
