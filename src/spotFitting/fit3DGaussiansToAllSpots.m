function Spots = fit3DGaussiansToAllSpots(Prefix, nSpots, varargin)
%%

cleanupObj = onCleanup(@myCleanupFun);

optionalResults = '';

segmentSpots = false;
displayFigures = false;
nWorkers = 8;
keepPool = false;
dogs = [];
saveType = '.tif';
save_flag = true;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'displayFigures')
        displayFigures = true;
    elseif strcmpi(varargin{i}, 'segmentSpots')
        Spots = varargin{i+1};
        segmentSpots = true;
    elseif strcmpi(varargin{i}, 'optionalResults')
        optionalResults = varargin{i+1};
    elseif strcmpi(varargin{i}, 'noSave')
        save_flag = false;
    elseif strcmpi(varargin{i}, 'nWorkers')
        nWorkers = varargin{i+1};
    elseif strcmpi(varargin{i}, 'keepPool')
        keepPool = true;
    elseif strcmpi(varargin{i}, 'dogs')
        dogs = varargin{i+1};
    elseif strcmpi(varargin{i}, 'saveAsMat') | strcmpi(varargin{i}, '.mat')
        saveType = '.mat';
    end
end


thisExperiment = liveExperiment(Prefix);

[~,ProcPath,DropboxFolder,~, PreProcPath,...
    ~, Prefix, ~,Channel1,Channel2,~, Channel3, spotChannels] = readMovieDatabase(Prefix, optionalResults);


DataFolder=[DropboxFolder,filesep,Prefix];

if ~segmentSpots
    Spots = getSpots(thisExperiment);
end

FrameInfo = getFrameInfo(thisExperiment);

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
%     parfor frame = 1:numFrames
    for frame = 1:numFrames %frames
        SpotsFr = SpotsCh(frame);

        nSpotsPerFrame = length(SpotsFr.Fits);
        for spot = 1:nSpotsPerFrame
            SpotsFr = fitSnip3D(SpotsFr, ch, spot, frame, thisExperiment, PreProcPath, FrameInfo, nSpots);
        end
        SpotsCh(frame) = SpotsFr;
        send(q, frame); %update the waitbar
    end
    
    if iscell(Spots) & length(Spots) > 1
        Spots{ch} = SpotsCh;
    else Spots = SpotsCh; end
    
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
    save([DataFolder,filesep,'Spots3DToken.mat'],'Spots3DToken', '-v6')
    disp('3D fitting done on all spots.')
    try close(waitbarFigure); end
    
end

    function nUpdateWaitbar(~)
        try waitbar(p/numFrames, waitbarFigure); end
        p = p + 1;
    end

end