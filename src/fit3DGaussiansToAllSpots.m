function Spots = fit3DGaussiansToAllSpots(prefix, varargin)
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

load([DataFolder,filesep,'FrameInfo.mat'], 'FrameInfo');

startParallelPool(nWorkers, displayFigures, keepPool);

if ~iscell(Spots)
    Spots = {Spots};
end


%%
for ch = spotChannels
    
    waitbarFigure = waitbar(0, ['Fitting 3D Gaussians: Channel ', num2str(ch)]);
    
    q = parallel.pool.DataQueue;
    afterEach(q, @nUpdateWaitbar);
    p = 1;
    
    SpotsCh = Spots{ch};
    numFrames = length(SpotsCh);
    
    
    for frame = 1:numFrames %frames
        nSpotsPerFrame = length(SpotsCh(frame).Fits);
        for spot = 1:nSpotsPerFrame
            SpotsCh = fitSnip3D(SpotsCh, ch, spot, frame, Prefix, PreProcPath, ProcPath, FrameInfo, dogs, displayFigures, saveType);
        end
        send(q, frame); %update the waitbar
        
    end
    
    Spots{ch} = SpotsCh;
    
end

if length(Spots) < 2
    Spots = Spots{1};
end


save([DataFolder,filesep,'Spots.mat'],'Spots', '-v7.3');
disp('3D fitting done on all spots.')

    function nUpdateWaitbar(~)
        waitbar(p/numFrames, waitbarFigure);
        p = p + 1;
    end

end