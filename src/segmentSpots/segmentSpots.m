% segmentSpots(Prefix, Threshold, [Options])
%
% DESCRIPTION
% Identify and segment individual transcription spots.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% Threshold: Threshold to be used (>=). For 8-bit images (color values go from 0 to 255)
%           should be kept at ~90-200 for lattice light-sheet data,
%           and at ~5-10 for confocal data (Leica SP8).
%           For 16-bit images (color values go from 0 to 65535) the numbers are orders
%           of magnitude larger, in the order of thousands.
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
% 'Shadows': This option should be followed by 0, 1 or 2. This
%            specifies the number of requisite z-planes in addition to the 
%            brightest plane that a spot must appear in for the spot to 
%            pass quality control. For shadows=1, that means each spot must
%            appear in at least 2 z-planes; for shadows=2, it's 3 z-planes
%            Default number of shadows is 1.
%
% 'keepPool': Don't shut down the parallel pool when the script is done
% running.
% 'nWorkers': Specify the number of workers to use during parallel
% processing
% 'noIntegralZ':  Don't establish center slice at position that maximizes raw fluo integral
%                 across sliding 3 z-slice window.
% 'autoThresh': Pops up a UI to help decide on a threshhold
% 'keepProcessedData': Keeps the ProcessedData folder for the given prefix after running segment spots
% 'fit3D': Fit 3D Gaussians to all segmented spots (assumes 1 locus per spot).
% 'fit3DOnly': Skip segmentation step and perform 3D fits
% 'skipChannel': Skips segmentation of channels inputted array (e.g. [1]
%                skips channel 1, [1, 2] skips channels 1 and 2
% 'optionalResults': use this if you have multiple Results/Dropbox folders
% for the same data to specify which you'll use.
% 'nuclearMask': Use the Ellipses structure to filter out particles
% detected outside of nuclei. 
%'track': track after running
% 'segmentChannel': use the DoGs of one channel to segment another
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



%% Argument validation


arguments   
    Prefix char
    Threshold (1,:) double
end

arguments (Repeating)
    varargin
end


[displayFigures, lastFrame, numShadows, keepPool, ...
    autoThresh, initialFrame, useIntegralCenter, Weka, keepProcessedData,...
    fit3D, skipChannel, optionalResults, filterMovieFlag, gpu, nWorkers, saveAsMat,...
    saveType, nuclearMask, DataType, track, skipSegmentation, frameRange, segmentChannel]...
    = determineSegmentSpotsOptions(varargin{:});

%validate the Threshold argument
if isempty(Threshold)
    Threshold = NaN;
end




%% Setup





cleanupObj = onCleanup(@myCleanupFun);
%this function uses persistent (static) variables to speed computation.
%if not cleared, this could lead to errors 
clear fitSingleGaussian
warning('off', 'MATLAB:MKDIR:DirectoryExists');




%% Main code




liveExperiment = LiveExperiment(Prefix);

spotChannels = liveExperiment.spotChannels;

[~, ~, DropboxFolder, ~, ~] = DetermineLocalFolders(Prefix, optionalResults);


if ~isempty(DataType)
     args = [Prefix, Threshold, varargin];
     writeScriptArgsToDataStatus(DropboxFolder, DataType, Prefix, args, 'Found filtered threshold', 'segmentSpots')
end

FrameInfo = getFrameInfo(liveExperiment);

DogOutputFolder=[liveExperiment.procFolder,filesep,'dogs',filesep];

if length(dir(DogOutputFolder)) <= 2
    error('Filtered movie files not found. Did you run FilterMovie?')
end

nSpotChannels = length(spotChannels);

%make sure the user inputted the right number of thresholds
if nSpotChannels > 1 && isempty(skipChannel) && length(Threshold) ~= nSpotChannels
    error('You must input the correct number of thresholds.');
end
   

if lastFrame==0
    lastFrame = numel(FrameInfo);
end
 
% The spot finding algorithm first segments the image into regions that are
% above the threshold. Then, it finds global maxima within these regions by searching in a region "neighborhood"
% within the regions.
pixelSize_nm = FrameInfo(1).PixelSize * 1000;
neighboorhood_nm = 1300;
neighborhood_px = round(neighboorhood_nm / pixelSize_nm);
snippetSize_px = 2 * (floor(neighboorhood_nm / (2 * pixelSize_nm))) + 1; % note that this is forced to be odd

falsePositives = 0;
if ~skipSegmentation
    disp('Segmenting spots...')
    Spots = cell(1, nSpotChannels);
    n = 0;
    for channelIndex = spotChannels
    
        n = n + 1;    
    
        if ismember(channelIndex, skipChannel)
            continue
        end
                
        if isempty(segmentChannel)
            segmentChannel = channelIndex;
        end

        tic;

        [ffim, doFF] = loadSegmentSpotsFlatField(...
            liveExperiment.userPreFolder, Prefix, spotChannels);
        
        tempSpots = segmentTranscriptionalLoci(channelIndex, initialFrame, lastFrame, FrameInfo(1).NumberSlices, ...
            liveExperiment.userPreFolder, Prefix, DogOutputFolder, displayFigures, doFF, ffim, Threshold(n), neighborhood_px, ...
            snippetSize_px, pixelSize_nm, FrameInfo(1).FileMode, [],...
             filterMovieFlag, optionalResults, gpu, saveAsMat, saveType, nuclearMask, autoThresh, segmentChannel);

        tempSpots = segmentSpotsZTracking(pixelSize_nm,tempSpots);

        [~, falsePositives, tempSpots] = findBrightestZ([], numShadows,...
            useIntegralCenter, 0, tempSpots);

        Spots{n} = tempSpots;

        timeElapsed = toc;
        disp(['Elapsed time: ', num2str(timeElapsed / 60), ' min'])
        try %#ok<TRYNC>
            log = logSegmentSpots(DropboxFolder, Prefix, timeElapsed, [], lastFrame, Spots, falsePositives, Threshold, channelIndex, numShadows, intScale, fit3D);
            display(log);
        end

    end
else
    disp('loading Spots.mat structure for 3D fitting...')
    load([DropboxFolder, filesep, Prefix, filesep, 'Spots.mat'], 'Spots');
end

%% Clean up Spots structure
% Remove duplicate spots 
% MT, 2022-06-08: not sure why identical duplicates are making their way
% into the Spots structure (with the TF cluster data), but it's best to fix
% this before we try to do 3D Gaussian fitting which is super time
% intensive
% [TODO]


%If we only have one channel, then convert Spots to a
%standard structure.
if nSpotChannels == 1 && iscell(Spots)
    Spots = Spots{1}; 
end

%% Save Spots.mat
mkdir([DropboxFolder, filesep, Prefix]);
if whos(var2str(Spots)).bytes < 2E9
    save([DropboxFolder, filesep, Prefix,...
        filesep, 'Spots.mat'], 'Spots', '-v6');
    disp('Spots.mat saved.')
else
    save([DropboxFolder, filesep, Prefix,...
        filesep, 'Spots.mat'], 'Spots', '-v7.3', '-nocompression');
    disp('Spots.mat saved.')
end

%% Track spots
if track, TrackmRNADynamics(Prefix, 'noretrack'); end

%% Fit 3D Gaussians to all segmented spots
if fit3D > 0
    disp('Fitting 3D Gaussians...')
    fit3DGaussiansToAllSpots(Prefix, 'segmentSpots', Spots, ...
                             'nWorkers', nWorkers);
    disp('3D Gaussian fitting completed.')
end

if ~keepPool    
    try  %#ok<TRYNC>
        poolobj = gcp('nocreate');
        delete(poolobj);
    end
end


end
