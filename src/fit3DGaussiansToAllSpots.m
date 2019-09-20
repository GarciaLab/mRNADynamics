function Spots = fit3DGaussiansToAllSpots(prefix, nSpots, varargin)
%%
optionalResults = '';

segmentSpots = false;
displayFigures = false;
nWorkers = 8;
keepPool = false;
dogs = [];
saveType = '.tif';

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'displayFigures')
        displayFigures = true;
    elseif strcmpi(varargin{i}, 'segmentSpots')
        Spots = varargin{i+1};
        segmentSpots = true;
    elseif strcmpi(varargin{i}, 'optionalResults')
        optionalResults = varargin{i+1};
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

[~,ProcPath,DropboxFolder,~, PreProcPath,...
    ~, Prefix, ~,Channel1,Channel2,~, Channel3, spotChannels] = readMovieDatabase(prefix, optionalResults);


DataFolder=[DropboxFolder,filesep,prefix];

if ~segmentSpots
    load([DataFolder,filesep,'Spots.mat'], 'Spots');
end

FrameInfo = load([DataFolder,filesep,'FrameInfo.mat'], 'FrameInfo');
FrameInfo = FrameInfo.FrameInfo;

% startParallelPool(nWorkers, displayFigures, keepPool);


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
    parfor frame = 1:numFrames %frames
        SpotsFr = SpotsCh(frame);

        nSpotsPerFrame = length(SpotsFr.Fits);
        for spot = 1:nSpotsPerFrame
            SpotsFr = fitSnip3D(SpotsFr, ch, spot, frame, Prefix, PreProcPath, FrameInfo, nSpots);
%             fitSnip3D(SpotsFr, spotChannel, spot, frame, Prefix, PreProcPath, FrameInfo)
        end
        SpotsCh(frame) = SpotsFr;
        send(q, frame); %update the waitbar
    end
    
    if iscell(Spots) & length(Spots) > 1
        Spots{ch} = SpotsCh;
    else
        Spots = SpotsCh;
    end
    
end

if iscell(Spots) & length(Spots) < 2
    Spots = Spots{1};
end


save([DataFolder,filesep,'Spots.mat'],'Spots', '-v7.3');
disp('3D fitting done on all spots.')
close(waitbarFigure);

    function nUpdateWaitbar(~)
        waitbar(p/numFrames, waitbarFigure);
        p = p + 1;
    end

end