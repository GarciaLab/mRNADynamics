function Spots = fit3DGaussiansToAllSpots(Prefix, varargin)
%%

cleanupObj = onCleanup(@myCleanupFun);

optionalResults = '';

segmentSpots = [];
displayFigures = false;
nWorkers = 8;
keepPool = false;
dogs = [];
save_flag = true;
nSpots = 1;


%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end
% 
% for i = 1:length(varargin)
%     if strcmpi(varargin{i}, 'displayFigures')
%         displayFigures = true;
%     elseif strcmpi(varargin{i}, 'segmentSpots')
%         Spots = varargin{i+1};
%         segmentSpots = true;
%     elseif strcmpi(varargin{i}, 'optionalResults')
%         optionalResults = varargin{i+1};
%     elseif strcmpi(varargin{i}, 'noSave')
%         save_flag = false;
%     elseif strcmpi(varargin{i}, 'nWorkers')
%         nWorkers = varargin{i+1};
%     elseif strcmpi(varargin{i}, 'keepPool')
%         keepPool = true;
%     elseif strcmpi(varargin{i}, 'dogs')
%         dogs = varargin{i+1};
%     end
% end

thisExperiment = liveExperiment(Prefix);

DataFolder = thisExperiment.resultsFolder;
PreProcPath = thisExperiment.preFolder;
spotChannels = thisExperiment.spotChannel;

if ~isempty(segmentSpots)
    Spots = getSpots(thisExperiment);
end

FrameInfo = getFrameInfo(thisExperiment);

movieMat = getMovieMat(thisExperiment);

startParallelPool(nWorkers, displayFigures, keepPool);


%%
for ch = spotChannels
    
    waitbarFigure = waitbar(0, ['Fitting 3D Gaussians: Channel ', num2str(ch)]);
    
    q = parallel.pool.DataQueue;
    afterEach(q, @nUpdateWaitbar);
    p = 1;
    
    if iscell(Spots)
        SpotsCh = Spots{ch};
    else
        SpotsCh = Spots;
    end
    
    numFrames = length(SpotsCh);
    
    % iterate through frames
    for frame = 1:numFrames %frames
        SpotsFr = SpotsCh(frame);
        
        nSpotsPerFrame = length(SpotsFr.Fits);
        for spotIndex = 1:nSpotsPerFrame
            SpotsFr = fitSnip3D(SpotsFr, ch, spotIndex, frame, Prefix, PreProcPath, FrameInfo, nSpots, movieMat);
        end
        SpotsCh(frame) = SpotsFr;
        send(q, frame); %update the waitbar
    end
    
    if iscell(Spots) && length(Spots) > 1
        Spots{ch} = SpotsCh;
    else
        Spots = SpotsCh;
    end
    
end

if iscell(Spots) & length(Spots) < 2
    Spots = Spots{1};
end

if save_flag
    
    if whos(var2str(Spots)).bytes < 2E9
        save([DataFolder,filesep,'Spots.mat'],'Spots', '-v6');
    else
        save([DataFolder,filesep,'Spots.mat'],'Spots', '-v7.3', '-nocompression');
    end
    
    Spots3DToken = now;
    save([DataFolder,filesep,'Spots3DToken.mat'],'Spots3DToken')
    disp('3D fitting done on all spots.')
    try close(waitbarFigure); end
    
end

    function nUpdateWaitbar(~)
        try waitbar(p/numFrames, waitbarFigure); end
        p = p + 1;
    end

end