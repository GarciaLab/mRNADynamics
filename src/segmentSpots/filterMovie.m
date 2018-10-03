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
% 'keepPool': Don't shut down the parallel pool when the script is done
% running.
% 'Frames', N:Run the code from frame 1 to frame N. Defaults to all
%                frames. It's suggested to run 5-20 frames for debugging.
% 'nWorkers': Specify the number of workers to use during parallel
% processing
%
% 'highPrecision': Uses higher precision filtering for segmentation when finding dogs
% 'customFilter': Choose which filter to use to segment the image. Name
%                  should be a string, followed by a cell with your filter
%                  or filters2018-09-20-pNos-FLP-7_embryo2_581power_2_0_Linear unmixing
%                 ex. filterMovie(Prefix,'customFilter', 'Structure_largest', {1, 8})
%           Filter Options:
%               'Gaussian_blur''Median'
%               'Maximum'
%               'Laplacian''Minimum'
%               'Mean'
%               'Hessian_largest''Hessian_smallest'
%               [DEFAULT] 'Difference_of_Gaussian' (2 sigmas) [DEFAULT]
%               
%
% OUTPUT
%
% Author (contact): Armando Reimer (areimer@berkeley.edu) / Matías Potel Feola (matias.potel.feola@gmail.com)
% Created: 8/23/2018
%
% Documented by: Armando Reimer (areimer@berkeley.edu)

function log = filterMovie(Prefix, varargin)

  warning('off', 'MATLAB:MKDIR:DirectoryExists');

  [displayFigures, numFrames, customFilter, highPrecision,...
      filterType, keepPool, sigmas, nWorkers, app, kernelSize] = determineFilterMovieOptions(varargin);

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
    else 
        try  %#ok<TRYNC>
            poolobj = gcp('nocreate');
            delete(poolobj);
        end
    end 
  
  clear rawdir;

  pixelSize = FrameInfo(1).PixelSize * 1000; %nm
  close all;

  coatChannel = getCoatChannel(ExperimentType, Channel1, Channel2);

  [sigmas] = generateDifferenceOfGaussianImages(DogOutputFolder, pixelSize, customFilter, nCh, ExperimentType, ...
    coatChannel, numFrames, displayFigures, zSize, PreProcPath, Prefix, filterType, highPrecision, sigmas, nWorkers, app,kernelSize);

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
  
  if ~keepPool

    try  %#ok<TRYNC>
      poolobj = gcp('nocreate');
      delete(poolobj);
    end 

  end 

end 
