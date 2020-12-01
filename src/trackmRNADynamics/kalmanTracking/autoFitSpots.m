function Spots = autoFitSpots(Prefix, autoFitInfo, Spots, varargin)
%%
cleanupObj = onCleanup(@myCleanupFun);

segmentSpots = false;
displayFigures = false;
nWorkers = 8;
keepPool = true;
% nSpots = 2;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'displayFigures')
        displayFigures = true;
    elseif strcmpi(varargin{i}, 'nWorkers')
        nWorkers = varargin{i+1};
    elseif strcmpi(varargin{i}, 'keepPool')
        keepPool = true;  
    end
end

disp(['auto-Fitting 3D Gaussians to: ', Prefix]);

% generate liveExperiment object
liveExperiment = LiveExperiment(Prefix);

% get mat containing raw image data
disp('Loading image files...')
movieMat = getMovieMat(liveExperiment);
disp('Done')

% create dircetory
DataFolder=[liveExperiment.resultsFolder,filesep];

% initialize pool
startParallelPool(nWorkers, displayFigures, keepPool);

for CurrentChannel = liveExperiment.spotChannels       
    
    if iscell(Spots)
        SpotsCh = Spots{CurrentChannel};
    else
        SpotsCh =  Spots;
    end
    
    if ~isempty(movieMat)
        movieMatCh = double(movieMat(:, :, :, :, CurrentChannel));
    else
        movieMatCh = [];
    end
    
    % get unique list of frames to fit
    FrameIndex = unique(autoFitInfo(CurrentChannel).Frame);
    
    Spots = autoFitting2D()
%     if ~isempty(movieMatCh)
%         parfor frame_i = 1:length(FrameIndex)
%             frame = FrameIndex(frame_i);
%             imStack = movieMatCh(:, :, :,frame);   
%             preSpots(frame_i) = spotFittingLoop(SpotsCh(frame).Fits(preIndexVec(preFrameVec==frame)), liveExperiment, imStack, [], nSpots);
%             % update waitbar
%             send(q, frame); %update the waitbar
%         end
%     else
%         parfor frame_i = 1:length(FrameIndex)
%             frame = FrameIndex(frame_i);
%             imStack = getMovieFrame(liveExperiment, frame, ch);
%             preSpots(frame_i) = spotFittingLoop(SpotsCh(frame).Fits(preIndexVec(preFrameVec==frame)), liveExperiment, imStack, [], nSpots);
%             % update waitbar
%             send(q, frame); %update the waitbar
%         end
%     end
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
