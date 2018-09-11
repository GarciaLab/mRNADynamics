% filterMovie(Prefix, [Options])
%
% DESCRIPTION
% Generates differences of gaussians files that will be later used to run segment spots.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% [Options]: See below.
%
% OPTIONS
% 'Frames', N:Run the code from frame 1 to frame N. Defaults to all
%                frames. It's suggested to run 5-20 frames for debugging.
%
% 'highPrecision': Uses higher precision filtering for segmentation when finding dogs
% 'customFilters': Choose which filter to use to segment the image. Name
%                  should be a string, followed by a cell with your filter
%                  or filters
%                 ex. filterMovie(Prefix,'customFilter', 'Structure_largest', {1, 8})
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
%
% Author (contact): Armando Reimer (areimer@berkeley.edu) / Matías Potel Feola (matias.potel.feola@gmail.com)
% Created: 8/23/2018
%
% Documented by: Armando Reimer (areimer@berkeley.edu)

function log = filterMovie(Prefix, varargin)

  warning('off', 'MATLAB:MKDIR:DirectoryExists');

  [displayFigures, numFrames, customFilter, highPrecision, filterType, sigmas] = determineFilterMovieOptions(varargin);

  % Start timer
  tic;

  [~, ~, ~, ~, ~, ~, ~, ExperimentType, Channel1, Channel2, ~] = readMovieDatabase(Prefix);

  [~, FISHPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

  load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat']);

  DogOutputFolder = [FISHPath, filesep, Prefix, '_', filesep, 'dogs'];
  mkdir(DogOutputFolder)

  zSize = FrameInfo(1).NumberSlices + 2;

  if numFrames == 0
    numFrames = length(FrameInfo);
  end 

  nCh = 1;

  if strcmpi(ExperimentType, '2spot2color')
    nCh = 2;
  end 

  clear rawdir;

  pixelSize = FrameInfo(1).PixelSize * 1000; %nm
  close all force;

  coatChannel = getCoatChannel(ExperimentType, Channel1, Channel2);

  [sigmas] = generateDifferenceOfGaussianImages(DogOutputFolder, pixelSize, customFilter, nCh, ExperimentType, ...
    coatChannel, numFrames, displayFigures, zSize, PreProcPath, Prefix, filterType, highPrecision, sigmas);

  t = toc;
  
  disp(['Elapsed time: ', num2str(t / 60), ' min'])

  logFile = [DropboxFolder, filesep, Prefix, filesep, 'log.mat'];

  if exist(logFile, 'file')
    load(logFile);
  else 
    log = struct();
  end 

  log(end + 1).Date = date;
  log(end).runTime = t / 60; % min
  log(end).LastFrame = numFrames;
  log(end).TimePerFrame = (t / 60) / numFrames;
  log(end).Filter = filterType;
  log(end).sigmas = sigmas;

  save(logFile, 'log', '-v7.3');

end 
