function Spots = fit3DGaussiansToAllSpots(Prefix, varargin)
%
% DESCRIPTION
%
%
% INPUT ARGUMENTS
% Prefix:
%
%
% OPTIONS
%
%
% OUTPUT
% Author (contact): Nicholas Lammers (nlammers@berkeley.edu)
% Created: 2021-03-04

cleanupObj = onCleanup(@myCleanupFun);

segmentSpots = false;
displayFigures = false;
nWorkers = 8;
keepPool = true;
save_flag = true;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'displayFigures')
        displayFigures = true;
    elseif strcmpi(varargin{i}, 'segmentSpots')
        
        %in the future, this should just be
        %replaced with a call to dbstack to check
        %if the caller is segmentSpots
       
        Spots = varargin{i+1};
        segmentSpots = true;  
    elseif strcmpi(varargin{i}, 'noSave')
        save_flag = false;
    elseif strcmpi(varargin{i}, 'nWorkers')
        nWorkers = varargin{i+1};
    elseif strcmpi(varargin{i}, 'keepPool')
        keepPool = true;  
%     elseif strcmpi(varargin{i}, 'fitSingleGaussian')
%         nSpots = 1;
    end
end

disp(['Fitting 3D Gaussians to: ', Prefix]);

% generate liveExperiment object
liveExperiment = LiveExperiment(Prefix);

% get mat containing raw image data
disp('Loading image files...')
movieMat = getMovieMat(liveExperiment);
disp('Done')

% create dircetory
DataFolder = [liveExperiment.resultsFolder,filesep];

if ~segmentSpots
    Spots = getSpots(liveExperiment);
end
startParallelPool(nWorkers, displayFigures, keepPool);

%% %%%%%%%%%%%%%%%%%%%%% Clean up Spots Structure %%%%%%%%%%%%%%%%%%%%%%%%%
% NL: implement this once I know which fields I'm adding
% error('pause here')

%% %%%%%%%%%%%%%%%%%%%%% Perform 3D fits to all spots %%%%%%%%%%%%%%%%%%%%%
for ch = liveExperiment.spotChannels       
    
    if iscell(Spots)
        SpotsCh = Spots{ch};
    else
        SpotsCh =  Spots;
    end
    
    if ~isempty(movieMat)
        movieMatCh = double(movieMat(:, :, :, :, ch));
    else
        movieMatCh = [];
    end
    
    numFrames = length(SpotsCh);
    
    % determine whether to expect 1 or 2 spots per locus
    % if this is TF cluster data, 1. For all transcription spots we assume
    % 2    
    inputChannels = liveExperiment.inputChannels;  
    nSpots = 2;
    if any(inputChannels==ch)
        nSpots = 1;
    end          
    
    waitbarFigure = waitbar(0, ['Fitting 3D Gaussians: Channel ', num2str(ch)]);
    
    q = parallel.pool.DataQueue;
    afterEach(q, @nUpdateWaitbar);
    p = 1;
    
    % iterate through all spots  
    if ~isempty(movieMatCh)
        parfor frame = 1:numFrames %frames                
            imStack = movieMatCh(:, :, :, frame);
            
            SpotsCh(frame) = spotFittingLoop(SpotsCh(frame).Fits, liveExperiment, imStack, nSpots);

            send(q, frame); %update the waitbar
        end
    else
        parfor frame = 1:numFrames %frames                
            imStack = getMovieFrame(liveExperiment, frame, ch);
            
            SpotsCh(frame) = spotFittingLoop(SpotsCh(frame).Fits, liveExperiment, imStack, nSpots);

            send(q, frame); %update the waitbar
        end
    end
    if iscell(Spots) && length(Spots) > 1
        Spots{ch} = SpotsCh;
    else
        Spots = SpotsCh; 
    end
     
end

if iscell(Spots) && length(Spots) < 2
    Spots = Spots{1};
end

if save_flag
    if whos(var2str(Spots)).bytes < 2E9
        save([DataFolder,filesep,'Spots.mat'],'Spots', '-v6');
    else
        save([DataFolder,filesep,'Spots.mat'],'Spots', '-v7.3', '-nocompression');
    end
    Spots3DToken = now;
    save([DataFolder,filesep,'Spots3DToken.mat'],'Spots3DToken', '-v6')
    disp('3D fitting done on all spots.')
    try close(waitbarFigure); end
    
end
    
function nUpdateWaitbarPre(~)
    try waitbar(p/numSamples, waitbarFigure); end
    p = p + 1;
end

function nUpdateWaitbar(~)
    try waitbar(p/numFrames, waitbarFigure); end
    p = p + 1;
end

end