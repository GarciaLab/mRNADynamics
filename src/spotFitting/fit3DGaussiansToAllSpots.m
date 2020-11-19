function Spots = fit3DGaussiansToAllSpots(Prefix, varargin)
%%

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
DataFolder=[liveExperiment.resultsFolder,filesep,Prefix];

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
    
    % iterate through spots and pull basic indexing info and fluorescence
    % stats
    frameRefVec = [];
    indexRefVec = [];
    spotFluoVec = [];
    for frame = 1:numFrames   
        for ind = 1:length(SpotsCh(frame).Fits)
            frameRefVec(end+1) = frame;
            indexRefVec(end+1) = ind;
            zIndex = SpotsCh(frame).Fits(ind).brightestZ==SpotsCh(frame).Fits(ind).z;
            spotFluoVec(end+1) = SpotsCh(frame).Fits(ind).FixedAreaIntensity(zIndex);
        end
    end
    
    % identify spots falling into the 3rd quartile
    q31 = prctile(spotFluoVec,50);
    q32 = prctile(spotFluoVec,75);
    q3Indices = find(spotFluoVec<=q32&spotFluoVec>q31);        
    
    % perform preliminary fitting to estimate spotPSF dims
    preSpots = struct('Fits',[]);
    spotDims = [];
    if length(q3Indices) >= 50
        % randomly select indices to use
        preIndices = randsample(q3Indices,min([250 length(q3Indices)]),false);
        % get unique list of frames
        preFrameVec = frameRefVec(preIndices);
        preIndexVec = indexRefVec(preIndices);
        preFrameIndex = unique(preFrameVec);
                
        waitbarFigure = waitbar(0, ['Performing initial fits to estimate PSF dimensions ', num2str(ch)]);
    
        q = parallel.pool.DataQueue;
        afterEach(q, @nUpdateWaitbarPre);
        p = 1;
        
        if ~isempty(movieMatCh)
            parfor frame = preFrameIndex(1:2)%numSamples
                imStack = movieMatCh(:, :, :,frame);   
                preSpots(frame) = spotFittingLoop(SpotsCh(frame).Fits(preIndexVec(preFrameVec==frame)), liveExperiment, imStack, []);
                % update waitbar
                send(q, frame); %update the waitbar
            end
        else
            parfor frame = preFrameIndex(1:10)%numSamples
                imStack = getMovieFrame(liveExperiment, frame, ch);
                preSpots(frame) = spotFittingLoop(SpotsCh(frame).Fits(preIndexVec(preFrameVec==frame)), liveExperiment, imStack, []);
                % update waitbar
                send(q, frame); %update the waitbar
            end
        end
        % extract parameters
        spotParamMat = [];
        preFits = [preSpots.Fits];
        for i = 1:length(preFits)                        
            spotParamMat = vertcat(spotParamMat, preFits(i).SpotFitInfo3D.RawFitParams);            
        end

        % calculate average inferred spot dimensions
        spotDims.sigmaXY = mean(spotParamMat(:,2));
        spotDims.sigmaXYSE = min([.5*spotDims.sigmaXY std(spotParamMat(:,2))]);
        spotDims.sigmaZ = mean(spotParamMat(:,3));
        spotDims.sigmaZSE = min([.5*spotDims.sigmaZ std(spotParamMat(:,3))]);
               
    end
    
    waitbarFigure = waitbar(0, ['Fitting 3D Gaussians: Channel ', num2str(ch)]);
    
    q = parallel.pool.DataQueue;
    afterEach(q, @nUpdateWaitbar);
    p = 1;
    
    % iterate through all spots    
    if ~isempty(movieMatCh)
        parfor frame = 1:numFrames %frames                
            imStack = movieMatCh(:, :, :, frame);
            
            SpotsCh(frame) = spotFittingLoop(SpotsCh(frame).Fits, liveExperiment, imStack, spotDims);

            send(q, frame); %update the waitbar
        end
    else
        parfor frame = 1:numFrames %frames                
            imStack = getMovieFrame(liveExperiment, frame, ch);
            
            SpotsCh(frame) = spotFittingLoop(SpotsCh(frame).Fits, liveExperiment, imStack, spotDims);

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