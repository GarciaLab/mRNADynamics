function [displayFigures, numFrames, numShadows, intScale, keepPool, threshGUI, initialFrame, useIntegralCenter, Weka, keepProcessedData, fit3D, skipChannel, optionalResults, filterMovieFlag] = determineSegmentSpotsOptions(varargin)

varargin = varargin{1};

% Default options
displayFigures = 0;
numFrames = 0;
numShadows = 2;
intScale = 1;
nWorkers = 8;
keepPool = 0;
threshGUI = 0;
useIntegralCenter = 1;
initialFrame = 1;
Weka = 0;
keepProcessedData = false;
fit3D = 0;
skipChannel = [];
optionalResults = '';
filterMovieFlag = false;


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
            numFrames = varargin{i + 1};
        end
        
    elseif strcmpi(varargin{i}, 'InitialFrame')

        if ~isnumeric(varargin{i + 1}) || varargin{i + 1} < 1
            error('Wrong input parameter for initial frame.')
        else
            initialFrame = varargin{i + 1};
        end
        
    elseif strcmpi(varargin{i}, 'keepPool')
        keepPool = 1;
    elseif strcmpi(varargin{i}, 'intScale')
        intScale = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'noIntegralZ')
        useIntegralCenter = 0;
    elseif strcmpi(varargin{i}, 'skipChannel')
        skipChannel = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'nWorkers')
        nWorkers = varargin{i + 1};        
    elseif strcmpi(varargin{i}, 'autoThresh')
        threshGUI = 1;
    elseif strcmpi(varargin{i}, 'fit3D')
        fit3D = 1;
    elseif strcmpi(varargin{i}, 'filterMovie')
        filterMovieFlag = true;
    elseif strcmpi(varargin{i}, 'Weka')
        Weka = 1;
    elseif strcmpi(varargin{i}, 'optionalResults')
        optionalResults = varargin{i+1};
    elseif strcmpi(varargin{i}, 'tifs')
        error('Tifs generation is no longer supported from segmentSpotsML, try filterMovie(Prefix, ''Tifs'') instead.');  
  
    elseif strcmpi(varargin{i}, 'keepProcessedData')
      keepProcessedData = true;  
    end
    
end

startParallelPool(nWorkers, displayFigures, keepPool);
    
end