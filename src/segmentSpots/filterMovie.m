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
%               'Gaussian_blur'
%               'Median'
%               'Maximum'
%               'Laplacian''Minimum'
%               'Mean'
%               'Hessian_largest''Hessian_smallest' 'bright_spot_psf'
%               [DEFAULT] 'Difference_of_Gaussian' (2 sigmas) [DEFAULT]
%
%               for the 3D versions of these filters, append '_3D' to the
%               end of the filtername. output is still 2D tifs.
% saveAsMat: save as a .mat file instead of .tif
% nogpu: do not use the gpu
%
%
% Author (contact): Armando Reimer (areimer@berkeley.edu) / Matías Potel Feola (matias.potel.feola@gmail.com)
% Created: 8/23/2018
% Last updated: 10/08/2018 - Matías Potel Feola - Refactor filterMovie and segmentSpotsML
%
% Documented by: Matías Potel Feola (harrypotel@gmail.com)
function log = filterMovie(Prefix, varargin)


processType = 'basic';


disp(['Generating filtered movie from ', Prefix,'...']);
warning('off', 'MATLAB:MKDIR:DirectoryExists');

cleanupObj = onCleanup(@myCleanupFun);

dogs = [];

% Start timer
tic;

liveExperiment = LiveExperiment(Prefix);

ExperimentType = liveExperiment.experimentType;
spotChannels = liveExperiment.spotChannels;
ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;
MS2CodePath = liveExperiment.MS2CodePath;
PreProcPath = liveExperiment.userPreFolder;

FrameInfo = getFrameInfo(liveExperiment);

[displayFigures, numFrames, initialFrame, highPrecision, filterType, keepPool,...
    sigmas, nWorkers, app, kernelSize, Weka, justTifs,...
    ignoreMemoryCheck, classifierFolder, ...
    classifierPathCh1, customML, noSave, numType, gpu,...
    saveAsMat, saveType, DataType] = determineFilterMovieOptions(FrameInfo,varargin);

if ~isempty(DataType)
    args = varargin;
    writeScriptArgsToDataStatus(DropboxFolder, DataType, Prefix, args, 'Made filtered spot channel files', 'filterMovie')
end

zSize = 2;
for i = 1:size(FrameInfo,2)
    if (FrameInfo(i).NumberSlices+2)>zSize
        zSize = FrameInfo(i).NumberSlices + 2;
    end
end

if numFrames == 0
    numFrames = length(FrameInfo);
end

nSpotChannels = length(spotChannels);

%generate tif stacks

stacksFolder = [PreProcPath, filesep, Prefix, filesep, 'stacks'];
stacksExist = exist(stacksFolder, 'dir') &&...
    ~isempty(dir([stacksFolder, filesep, '*.tif']));

if (Weka || justTifs) && ~stacksExist
    generateTifsForWeka(Prefix, PreProcPath, numFrames,...
        nSpotChannels,spotChannels, zSize, initialFrame);
end

if Weka
    processType = 'weka';
elseif customML
    processType = 'customML';
end

if ~justTifs
    
    switch processType
        
        case  'basic'
            
            generateDifferenceOfGaussianImages(ProcPath,...
                spotChannels,...
                numFrames, displayFigures, zSize, PreProcPath,...
                Prefix, filterType, highPrecision, sigmas, app,...
                kernelSize, noSave, numType, gpu, saveAsMat, saveType);
            
        case 'weka'
            
            generateDogsWeka(Prefix, ProcPath, MS2CodePath,...
                PreProcPath, spotChannels, zSize, numFrames, nSpotChannels,...
                initialFrame, ignoreMemoryCheck, classifierPathCh1, classifierFolder);
            
        case 'customML'
            
            generateProbMapsCustomML(Prefix, ProcPath,...
                MS2CodePath, PreProcPath, ExperimentType, coatChannel, zSize, numFrames, nSpotChannels,...
                initialFrame, ignoreMemoryCheck, classifierPathCh1, classifierFolder);
            
        otherwise
            
            error('Processing type not recognized.')
            
    end
    
end

t = toc;

disp(['Elapsed time: ', num2str(t / 60), ' min'])

if ~justTifs
    try log = writeFilterMovieLog(t, Weka, DropboxFolder, Prefix,...
            initialFrame, numFrames, filterType, sigmas, classifierPathCh1);
    end
end

if ~keepPool && ~Weka && ~justTifs
    try %#ok<TRYNC>
        poolobj = gcp('nocreate');
        delete(poolobj);
    end
end

disp([Prefix, ' filtered.']);

end
