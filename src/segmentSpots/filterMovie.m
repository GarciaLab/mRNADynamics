% filterMovie(Prefix, [Options])
%
% DESCRIPTION
% Generates differences of gaussians files that will be later used to run segment spots using regular method or Weka.
% Also generates Tifs to create Weka classifier.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% [Options]: See below.
%
% OPTIONS
% 'Weka': Generates differences of gaussians files that will be later used to run segment spots ML, 
%         requires a classifier.
%
% 'Tifs' : Generates the TIF stacks necessary for doing Weka classification.
%          Recommended to run this before making a new classifier.
%
% 'displayFigures':   If you want to display plots and images.
%
% 'Frames', N: Run the code from frame 1 to frame N. Defaults to all
%              frames. It's suggested to run 5-20 frames for debugging.
%
% 'InitialFrame', N: Run the code from frame N to last frame. Defaults to first
%                    frame.
%                 Only for Weka and Tifs modes.
% 
% 'LastFrame', M: Run the code from initial frame to frame M. Defaults to all
%                 frames. It's suggested to run 5-20 frames for debugging.
%                 Only for Weka and Tifs modes.
%
% 'keepPool': Don't shut down the parallel pool when the script is done
%             running.
%
%
% 'nWorkers': Specify the number of workers to use during parallel
%             processing
%
% 'highPrecision': Uses higher precision filtering for segmentation when finding dogs
%
% 'ignoreMemoryCheck' : Ignores the memory check, should be used only for testing. Only for Weka.
% 'customML' : run Nick's custom ML
%
% 'customFilter': Choose which filter to use to segment the image. Name
%                  should be a string, followed by a cell with your filter
%                  or filters
%                 ex. filterMovie(Prefix,'customFilter', 'Structure_largest', {1, 8})
%           Filter Options:
%               'Gaussian_blur''Median'
%               'Maximum'
%               'Laplacian''Minimum'
%               'Mean'
%               'Hessian_largest''Hessian_smallest' 'bright_spot_psf'
%               [DEFAULT] 'Difference_of_Gaussian' (2 sigmas) [DEFAULT]
%
%               for the 3D versions of these filters, append '_3D' to the
%               end of the filtername. output is still 2D tifs. 
%
%
% Author (contact): Armando Reimer (areimer@berkeley.edu) / Matías Potel Feola (matias.potel.feola@gmail.com)
% Created: 8/23/2018
% Last updated: 10/08/2018 - Matías Potel Feola - Refactor filterMovie and segmentSpotsML
%
% Documented by: Matías Potel Feola (harrypotel@gmail.com)
function log = filterMovie(Prefix, varargin)

  warning('off', 'MATLAB:MKDIR:DirectoryExists');

  [displayFigures, numFrames, initialFrame, customFilter, highPrecision, filterType, keepPool,...
    sigmas, nWorkers, app, kernelSize, Weka, justTifs, ignoreMemoryCheck, classifierFolder, ...
    classifierPathCh1, customML] = determineFilterMovieOptions(varargin);

  % Start timer
  tic;
  
  startParallelPool(nWorkers, displayFigures);

  [~, ~, ~, ~, ~, ~, ~, ExperimentType, Channel1, Channel2, ~] = readMovieDatabase(Prefix);

  [~, ProcPath, DropboxFolder, MS2CodePath, PreProcPath] = DetermineLocalFolders(Prefix);

  
  load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat'], 'FrameInfo');
  zSize = FrameInfo(1).NumberSlices + 2;

  if numFrames == 0
    numFrames = length(FrameInfo);
  end

  nCh = 1;

  if strcmpi(ExperimentType, '2spot2color')
    nCh = 2;
  end

  coatChannel = getCoatChannel(ExperimentType, Channel1, Channel2);
  
  if ~Weka && ~justTifs
    sigmas = generateDifferenceOfGaussianImages(ProcPath, customFilter, nCh, ExperimentType, FrameInfo, coatChannel,...
      numFrames, displayFigures, zSize, PreProcPath, Prefix, filterType, highPrecision, sigmas, nWorkers, app, kernelSize);
  elseif Weka
      if ~exist([PreProcPath, filesep, Prefix, filesep, 'stacks'], 'dir')
        generateTifsForWeka(Prefix, ExperimentType, PreProcPath, numFrames, nCh,coatChannel, zSize, initialFrame);  
      end
    generateDogsWeka(Prefix, ProcPath, MS2CodePath, PreProcPath, ExperimentType, coatChannel, zSize, numFrames, nCh,...
      initialFrame, ignoreMemoryCheck, classifierPathCh1, classifierFolder);
  elseif justTifs
    generateTifsForWeka(Prefix, ExperimentType, PreProcPath, numFrames, nCh,coatChannel, zSize, initialFrame);
  elseif customML
      generateProbMapsCustomML(Prefix, ProcPath, MS2CodePath, PreProcPath, ExperimentType, coatChannel, zSize, numFrames, nCh,...
      initialFrame, ignoreMemoryCheck, classifierPathCh1, classifierFolder);
  end

  t = toc;

  disp(['Elapsed time: ', num2str(t / 60), ' min'])

  if ~justTifs
    log = writeFilterMovieLog(t, Weka, DropboxFolder, Prefix, initialFrame, numFrames, filterType, sigmas, classifierPathCh1);
  end
  
  if ~keepPool && ~Weka && ~justTifs
    try %#ok<TRYNC>
      poolobj = gcp('nocreate');
      delete(poolobj);
    end
  end

end
