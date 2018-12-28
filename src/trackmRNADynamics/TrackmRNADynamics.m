%[Particles,schnitzcells]=TrackmRNADynamics(varargin)
%
% DESCRIPTION
% %This function tracks transcription loci over time after
% segmentation and z-tracking have been performed. If nuclei have been
% tracked it uses that information for the particle tracking.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% Threshold1 : Primary. This must be an array of two thresholds for 2 spot 2
% color experiments.
% Threshold2 : Secondary. This must be an array of two thresholds for 2 spot 2
% color experiments.
% [Options]: See below.
%
% OPTIONS
% None.
%
% OUTPUT
% Particles.mat : List of time traces found in the movie.
% *_lin.mat : Schnitzcell nuclear tracking is modified here.
%
% Author (contact): Hernan Garcia (hggarcia@berkeley.edu)
% Created: 01/01/2013 ish.
% Last Updated: 9/11/2018
%
% Documented by: Armando Reimer (areimer@berkeley.edu)
%
%
%
% To do: The no-histone part of the code doesn't take into account the
% Approved field of the Particles structure.
% ^ AR 9/3/18: has this been done?
%
function [Particles, schnitzcells] = TrackmRNADynamics(varargin)

  [~, ~, DefaultDropboxFolder, ~, ~] = DetermineLocalFolders;

  [Prefix, Threshold1, Threshold2, app, bypassUserPrompt] = parseTrackmRNADynamicsArguments(DefaultDropboxFolder, varargin{:});

  % Save a backup of the threshold. We'll check that it didn't change if we're
  % doing retracking. This is definitely not an elegant solution
  Threshold1Backup = Threshold1;
  Threshold2Backup = Threshold2;

  % Get the actual folder now that we have the Prefix
  [~, ~, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

  % What type of experiment are we dealing with? Get this out of MovieDatabase
  [~, ExperimentType, ~, ~, ~, ~, Channel1, Channel2, ~, ~, ~, ~, ~, ...
      nc9, nc10, nc11, nc12, nc13, nc14, ~] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

  % Set the destination folder
  OutputFolder = [DropboxFolder, filesep, Prefix];

  % Determine the search radius based on the imaging conditions
  SearchRadiusMicrons = 3; % Search radius in um

  [Threshold1, Threshold2] = deleteTrackingIfThresholdsChangedRetracking(Threshold1, Threshold2, Threshold1Backup, ...
    Threshold2Backup, OutputFolder, bypassUserPrompt);

  % Load the information about this image
  % Check if we have FrameInfo otherwise try to get the information straight
  % from the file.
  [FrameInfo, PixelSize] = obtainPixelSize(OutputFolder);

  SearchRadius = ceil(SearchRadiusMicrons / PixelSize);

  % Check if we have tracked the lineages of the nuclei
  if exist([DropboxFolder, filesep, Prefix, filesep, Prefix, '_lin.mat'], 'file')
    UseHistone = 1;

    % Load the nuclei segmentation information
    load([DropboxFolder, filesep, Prefix, filesep, 'Ellipses.mat'], 'Ellipses')

    % Load the nuclei tracking information
    load([DropboxFolder, filesep, Prefix, filesep, Prefix, '_lin.mat'])

    % Do a bunch of check on schnitzcells only if tracking for the first
    % time.

  else
    UseHistone = 0;
    warning('Warning: No nuclei lineage tracking found. Proceeding with tracking particles only.')
  end

  validateExperimentTypeSupported(ExperimentType);

  SpotsChannel = determineSpotsChannel(ExperimentType, Channel1, Channel2);
  % Get the number of channels
  NCh = length(SpotsChannel);

  Particles = loadParticlesAndSelectForRetracking(OutputFolder, NCh);

  [Spots, SpotFilter] = loadSpotsAndCreateSpotFilter(DropboxFolder, Prefix, NCh, Threshold1);

  [ParticlesFig, particlesAxes, NucleiFig, nucAxes] = numberFoundParticles(app, UseHistone);
  NDigits = adjustIndexSizeAccordingToFrames(FrameInfo);

  [Particles, SpotFilter] = performTracking(Particles, schnitzcells, NCh, Spots, app, SpotFilter, PreProcPath, ...
    Prefix, UseHistone, ParticlesFig, SpotsChannel, NDigits, NucleiFig, particlesAxes, nucAxes, Ellipses, ...
    PixelSize, SearchRadius, ExperimentType, FrameInfo);
  
  mkdir([OutputFolder, filesep]);

  save([OutputFolder, filesep, 'Particles.mat'], 'Particles', 'SpotFilter', ...
    'Threshold1', 'Threshold2', '-v7.3');

  createFieldNCAndSaveFrameInfo(FrameInfo, OutputFolder, nc9, nc10, nc11, nc12, nc13, nc14);
end

% Check whether we're using the same threshold in the case of retracking
function [Threshold1, Threshold2] = deleteTrackingIfThresholdsChangedRetracking(Threshold1, Threshold2, Threshold1Backup, Threshold2Backup, OutputFolder, bypassUserPrompt)
  ParticlesFilePath = [OutputFolder, filesep, 'Particles.mat'];

  if exist(ParticlesFilePath, 'file')
    load(ParticlesFilePath, 'Threshold1', 'Threshold2')

    if (~sum(Threshold1 == Threshold1Backup) == length(Threshold1)) & ...
        (~sum(Threshold2 == Threshold2Backup) == length(Threshold2))

      if ~bypassUserPrompt
        Answer = input('Thresholds changed, will delete previous tracking. Proceed? (y/n):', 's');
      else
        Answer = 'y';
      end

      if strcmpi(Answer, 'y')
        Threshold1 = Threshold1Backup;
        Threshold2 = Threshold2Backup;
        delete(ParticlesFilePath)
      else
        error('Cannot retrack if the threshold changed')
      end

    end

  end

end

function [FrameInfo, PixelSize] = obtainPixelSize(OutputFolder)
  FrameInfoPath = [OutputFolder, filesep, 'FrameInfo.mat'];

  if exist(FrameInfoPath, 'file')
    load(FrameInfoPath, 'FrameInfo');

    %See if this came from the 2-photon, which is the default
    if ~isfield(FrameInfo, 'FileMode') || strcmp(FrameInfo(end).FileMode, 'TIF')
      PixelSize = 0.22;
    elseif strcmp(FrameInfo(1).FileMode, 'LSM') | strcmp(FrameInfo(1).FileMode, 'LSMExport')%#ok<*OR2>
      PixelSize = FrameInfo(1).PixelSize;
    elseif strcmp(FrameInfo(1).FileMode, 'LIFExport') || strcmp(FrameInfo(1).FileMode, 'LAT') || strcmp(FrameInfo(1).FileMode, 'DSPIN')%CS20170907
      PixelSize = FrameInfo(1).PixelSize;
    end

  else
    warning('No FrameInfo.mat detected. Trying to pull out magnification information from the TIF file')

    DZoom = dir([PreProcPath, filesep, Prefix, filesep, '*z*.tif']);
    ImageInfo = imfinfo([PreProcPath, filesep, Prefix, filesep, DZoom(1).name]);
    PixelSize = 1 / ImageInfo.XResolution;
  end

end

function validateExperimentTypeSupported(ExperimentType)

  if ~(strcmpi(ExperimentType, '1spot') || strcmpi(ExperimentType, '2spot') || ...
      strcmpi(ExperimentType, 'inputoutput') || ...
      strcmpi(ExperimentType, '2spot2color') || ...
      strcmpi(ExperimentType, 'lattice'))

    error('Experiment type in MovieDatabase not recognized')

  end

end

function SpotsChannel = determineSpotsChannel(ExperimentType, Channel1, Channel2)

  if strcmp(ExperimentType, '1spot') || strcmp(ExperimentType, '2spot')
    SpotsChannel = 1;
    %(MT, 2018-02-11) Added support for lattice imaging, maybe temporary -
    %FIX LATER
  elseif strcmp(ExperimentType, 'inputoutput') || strcmp(ExperimentType, 'lattice')
    SpotsChannel = find(~cellfun(@isempty, strfind(lower([Channel1, Channel2]), 'mcp')) | ...
      ~cellfun(@isempty, strfind(lower([Channel1, Channel2]), 'pcp')));

    if length(SpotsChannel) > 1
      error('only one output channel currently supported in inputoutput mode')
    end

  elseif strcmpi(ExperimentType, '2spot2color')
    SpotsChannel = [1, 2];
  end

end

function Particles = loadParticlesAndSelectForRetracking(OutputFolder, NCh);
  % Check if particle tracking has already been done on this dataset
  if exist([OutputFolder, filesep, 'Particles.mat'], 'file')

    load([OutputFolder, filesep, 'Particles.mat'], 'Particles')

    % If there's only one channel, Particles, Spots and other structures are
    % not saved as cells. We turn them into a cell to ensure
    % compatibility.
    if NCh == 1
      Particles = {Particles};
    end

    for Channel = 1:NCh

      if isfield(Particles{1}, 'Approved')
        display(['Performing retracking on channel ', num2str(Channel)])

        %Only keep the approved particles and start the tracking from there
        k = 1;

        for i = 1:length(Particles{Channel})

          if Particles{Channel}(i).Approved ~= 0
            NewParticles{Channel}(k) = Particles{Channel}(i);
            k = k + 1;
          end

        end

        if exist('NewParticles')
          Particles{Channel} = NewParticles{Channel};
        else
          Particles{Channel} = [];
        end

      else
        Particles{Channel} = [];
      end

    end

  else

    for Channel = 1:NCh
      Particles{Channel} = []; % This is the structure where we'll be tracking all particles.
    end

  end

end

function [Spots, SpotFilter] = loadSpotsAndCreateSpotFilter(DropboxFolder, Prefix, NCh, Threshold1)
  % First, generate a structure array with a flag that determines whether
  % each spots is above Threshold1
  if ~exist('Spots', 'var')
    load([DropboxFolder, filesep, Prefix, filesep, 'Spots.mat'], 'Spots')

    % If there's only one channel, Particles, Spots and other structures are
    % not saved as cells. We turn them into a cell to ensure
    % compatibility.
    if NCh == 1
      Spots = {Spots};
    end

    MaxSpots = cell(NCh);

    for Channel = 1:NCh
      %Determine the maximum number of spots in a given frame for the
      %whole movie
      MaxSpots{Channel} = 0;

      for i = 1:length(Spots{Channel})
        MaxSpots{Channel} = max([MaxSpots{Channel}, length(Spots{Channel}(i).Fits)]);
      end

      %This filter tells us whether a spot is above the threshold.
      if exist('SpotFilter', 'var')

        if ~iscell(SpotFilter)
          SpotFilter = {SpotFilter};
        end

      end

      SpotFilter{Channel} = nan(length(Spots{Channel}), MaxSpots{Channel});
      %Populate the filter
      for i = 1:length(Spots{Channel})

        for j = 1:length(Spots{Channel}(i).Fits)

          if sum(Spots{Channel}(i).Fits(j).DOGIntensity > Threshold1(Channel))
            SpotFilter{Channel}(i, j) = 1;
          end

        end

      end

    end

  end

end

function [ParticlesFig, particlesAxes, NucleiFig, nucAxes] = numberFoundParticles(app, UseHistone)
  %Start by numbering the particles found
  if isempty(app)
    ParticlesFig = figure;
    particlesAxes = axes(ParticlesFig);

    if UseHistone
      NucleiFig = figure;
      set(NucleiFig, 'units', 'normalized', 'position', [0.65, .5, .2, .2])
      nucAxes = axes(NucleiFig);
    end

  end
end

% See how  many frames we have and adjust the index size of the files to
% load accordingly
function NDigits = adjustIndexSizeAccordingToFrames(FrameInfo)

  if length(FrameInfo) < 1E3
    NDigits = 3;
  elseif length(FrameInfo) < 1E4
    NDigits = 4;
  else
    error('No more than 10,000 frames supported. Change this in the code')
  end

end

function createFieldNCAndSaveFrameInfo(FrameInfo, OutputFolder, nc9, nc10, nc11, nc12, nc13, nc14)
  % creating the field nc for FrameInfo
  if exist([OutputFolder, filesep, 'FrameInfo.mat'], 'file')
    numberOfFrames = length(FrameInfo);

    for currentFrame = 1:numberOfFrames

      if currentFrame < nc9
        FrameInfo(currentFrame).nc = 8;
      elseif (currentFrame >= nc9) & (currentFrame < nc10)
        FrameInfo(currentFrame).nc = 9;
      elseif (currentFrame >= nc10) & (currentFrame < nc11)
        FrameInfo(currentFrame).nc = 10;
      elseif (currentFrame >= nc11) & (currentFrame <= nc12)
        FrameInfo(currentFrame).nc = 11;
      elseif (currentFrame >= nc12) & (currentFrame <= nc13)
        FrameInfo(currentFrame).nc = 12;
      elseif (currentFrame >= nc13) & (currentFrame <= nc14)
        FrameInfo(currentFrame).nc = 13;
      elseif currentFrame >= nc14
        FrameInfo(currentFrame).nc = 14;
      end

    end

    save([OutputFolder, filesep, 'FrameInfo.mat'], 'FrameInfo')
  else
    warning('Tried to save nc frame information, but could not since there is no FrameInfo.mat')
  end
end