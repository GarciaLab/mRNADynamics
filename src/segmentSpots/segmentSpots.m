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
% 'autoThresh': Pops up a UI to help decide on a threshhold
% 'keepProcessedData': Keeps the ProcessedData folder for the given prefix after running segment spots
% 'fit3D': Fit 3D Gaussians to all segmented spots (assumes 1 locus per spot).
% 'fit3D2Spot': Fit 3D Gaussians to all segmented spots (assumes 2 loci per spot).
% 'skipChannel': Skips segmentation of channels inputted array (e.g. [1]
%                skips channel 1, [1, 2] skips channels 1 and 2
% 'optionalResults': use this if you have multiple Results/Dropbox folders
% for the same data to specify which you'll use.
% 'nuclearMask': Use the Ellipses structure to filter out particles
% detected outside of nuclei. 
%'track': track after running
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

[displayFigures, numFrames, numShadows, keepPool, ...
    autoThresh, initialFrame, useIntegralCenter, Weka, keepProcessedData,...
    fit3D, skipChannel, optionalResults, filterMovieFlag, gpu, nWorkers, saveAsMat,...
    saveType, nuclearMask, DataType, track]...
    = determineSegmentSpotsOptions(varargin{:});

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

[~, ~, ~, ~, ~, ~, ~, ~, ~, ~,~, ~, spotChannels] = readMovieDatabase(Prefix);

[~, ProcPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix, optionalResults);

if ~isempty(DataType)
     args = [Prefix, Threshold, varargin];
     writeScriptArgsToDataStatus(DropboxFolder, DataType, Prefix, args, 'Found filtered threshold', 'segmentSpots')
end

load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat'], 'FrameInfo');
if nuclearMask
    load([DropboxFolder, filesep, Prefix, filesep, 'Ellipses.mat'], 'Ellipses');
else
    Ellipses = {};
end


ProcessedDataFolder = [ProcPath, filesep, Prefix, '_'];
DogOutputFolder = [ProcessedDataFolder, filesep, 'dogs'];

microscope = FrameInfo(1).FileMode;

zSize = FrameInfo(1).NumberSlices;

nCh = length(spotChannels);

m = matfile([ProcPath, filesep, Prefix,  '_', filesep, Prefix, '_dogMat.mat']);
numFrames = size(m, 'dogMat', 1);

% The spot finding algorithm first segments the image into regions that are
% above the threshold. Then, it finds global maxima within these regions by searching in a region "neighborhood"
% within the regions.

pixelSize = FrameInfo(1).PixelSize * 1000; %nm
neighboorhood_size = 1300;
neighborhood = round(neighboorhood_size / pixelSize); %nm
snippet_size = 2 * (floor(neighboorhood_size / (2 * pixelSize))) + 1; % nm. note that this is forced to be odd
coatChannel = spotChannels;

falsePositives = 0;
Spots = cell(1, nCh);

for channelIndex = 1:nCh
    if ismember(channelIndex, skipChannel)
        continue
    end
    
    tic;
    
    [ffim, doFF] = loadSegmentSpotsFlatField(PreProcPath, Prefix, spotChannels); %let's not
    doFF = false;

    [tempSpots, dogs] = segmentTranscriptionalLoci(nCh, coatChannel, channelIndex, initialFrame, numFrames, zSize, ...
        PreProcPath, Prefix, ProcessedDataFolder, displayFigures, doFF, ffim, Threshold(channelIndex), neighborhood, ...
        snippet_size, pixelSize, microscope, Weka,...
         filterMovieFlag, optionalResults, gpu, saveAsMat, saveType, Ellipses);

    tempSpots = segmentSpotsZTracking(pixelSize,tempSpots);

    [~, falsePositives, tempSpots] = findBrightestZ([], numShadows, useIntegralCenter, 0, tempSpots, 'dogs', dogs);
                        
    Spots{channelIndex} = tempSpots;
    
    timeElapsed = toc;
    disp(['Elapsed time: ', num2str(timeElapsed / 60), ' min'])
    try %#ok<TRYNC>
        log = logSegmentSpots(DropboxFolder, Prefix, timeElapsed, [], numFrames, Spots, falsePositives, Threshold, channelIndex, numShadows, intScale, fit3D);
        display(log);
    end
    
end

%If we only have one channel, then convert Spots to a
%standard structure.
if nCh == 1 && iscell(Spots)
    Spots = Spots{1}; 
end

mkdir([DropboxFolder, filesep, Prefix]);
save([DropboxFolder, filesep, Prefix, filesep, 'Spots.mat'], 'Spots', '-v7.3');

if fit3D > 0
    disp('Fitting 3D Gaussians...')
    fit3DGaussiansToAllSpots(Prefix, fit3D, 'segmentSpots', Spots, 'nWorkers', nWorkers, saveType);
    disp('3D Gaussian fitting completed.')
end

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

if track
    TrackmRNADynamics(Prefix, 'noretrack');
end

end
