function [displayFigures, lastFrame, numShadows, keepPool, threshGUI, initialFrame, ...
    useIntegralCenter, Weka, keepProcessedData, fit3D, skipChannel,...
    optionalResults, filterMovieFlag, gpu, nWorkers, saveAsMat, saveType, ...
    nuclearMask, dataType, runTrackmRNADynamics, skipSegmentation, frameRange]...
    ...
    = determineSegmentSpotsOptions(varargin)

% Default options
displayFigures = false;
lastFrame = 0;
numShadows = 1;
nWorkers = 8;
keepPool = true;
threshGUI = false;
useIntegralCenter = true;
initialFrame = 1;
Weka = false;
keepProcessedData = true;
fit3D = 1;
skipSegmentation = false;
skipChannel = [];
optionalResults = '';
filterMovieFlag = false;
gpu = '';
saveAsMat = true;
saveType = '.mat';
nuclearMask = false;
dataType = '';
runTrackmRNADynamics = true;

  
for i = 1:length(varargin)
    
    if strcmpi(varargin{i}, 'displayFigures')
        displayFigures = 1;
        close all;
    elseif strcmpi(varargin{i}, 'Shadows')
        
        if (i + 1) > length(varargin) || ~ isnumeric(varargin{i + 1}) || varargin{i + 1} > 2
            error('Wrong input parameters. After ''Shadows'' you should input number of shadows(0, 1 or 2)')
        else
            numShadows = varargin{i + 1};
        end
        
    elseif strcmp(varargin{i}, 'Frames') || strcmpi(varargin{i}, 'LastFrame')
        
        if ~isnumeric(varargin{i + 1}) || varargin{i + 1} < 1
            error('Wrong input parameters. After ''Frames'' you should input the number of frames')
        else
            lastFrame = varargin{i + 1};
        end
        
    elseif strcmpi(varargin{i}, 'InitialFrame')

        if ~isnumeric(varargin{i + 1}) || varargin{i + 1} < 1
            error('Wrong input parameter for initial frame.')
        else
            initialFrame = varargin{i + 1};
        end
        
    elseif strcmpi(varargin{i}, 'keepPool')
        keepPool = true;
     elseif strcmpi(varargin{i}, 'dataSet') || strcmpi(varargin{i}, 'dataType')
        dataType = varargin{i+1};
    elseif strcmpi(varargin{i}, 'saveAsMat') || strcmpi(varargin{i}, '.mat')
        saveAsMat = true;
        saveType = '.mat';
     elseif strcmpi(varargin{i}, 'noGPU')
       gpu = 'noGPU';
     elseif strcmpi(varargin{i}, 'track')
       runTrackmRNADynamics = true;
    elseif strcmpi(varargin{i}, 'noIntegralZ')
        useIntegralCenter = 0;
    elseif strcmpi(varargin{i}, 'skipChannel')
        skipChannel = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'nWorkers')
        nWorkers = varargin{i + 1};        
    elseif strcmpi(varargin{i}, 'nuclearMask')
        error('.\n%s', ...
              'Error: The ''nuclearMask'' option is currently not compatible with the parfor loop in segmentTranscriptionalLoci.m (called by segmentSpots) due to an elusive, currently unsolved bug',...
              'To use this option, do the following:',...
              '(1) In segmentTranscriptionalLoci.m, change the parfor loop to a for loop (just delete the letters ''par'') at approx. line 127',...
              '(2) In segmentTranscriptionalLoci.m, find the two ''if shouldMaskNuclei'' statements, at approx. lines 95 and 198. Uncomment both these if statements in their entirety',...
              '(3) In determineSegmentSpotsOptions.m, comment out this error message at approx. line 82',...
              '(4) DO NOT push the changes in 1-3 to the master branch. If you need them, please stash on a separate branch',...
              'If you are confused or need help implementing the steps above, please contact Meghan Turner (meghan_turner@berkeley.edu).');
        if islogical(varargin{i+1})
            nuclearMask = varargin{i+1};
        else
            nuclearMask = true;
        end
    elseif strcmpi(varargin{i}, 'autoThresh')...
        || strcmpi(varargin{i}, 'determineThreshold')
        threshGUI = 1;
    elseif strcmpi(varargin{i}, 'fit3D')
        fit3D = 1;
    elseif strcmpi(varargin{i}, 'fit3DOnly')
        fit3D = 1;
        skipSegmentation = 1;
    elseif strcmpi(varargin{i}, 'fit3D2Spot') && fit3D == 0
        fit3D = 2;
    elseif strcmpi(varargin{i}, 'filterMovie')
        filterMovieFlag = true;
    elseif strcmpi(varargin{i}, 'optionalResults')
        optionalResults = varargin{i+1};
    elseif strcmpi(varargin{i}, 'keepProcessedData')
      keepProcessedData = true;  
    end
    
end

frameRange = [initialFrame, lastFrame]; 

startParallelPool(nWorkers, displayFigures, keepPool);
    
end