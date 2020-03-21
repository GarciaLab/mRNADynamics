function [movieMat, hisMat, maxMat, medMat, midMat]...
    = makeMovieMats(Prefix, PreProcPath, nWorkers, FrameInfo, varargin)

makeHis = false;
makeMovie = false;
makeProjs = false;

loadHis = true;
loadMovie = true;
loadProjs = true;

nPadding = 2; %normally we pad a blank image above and below the stack.
%if a movie is unpadded this will require modification


%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

%don't load things that aren't requested as outputs
if nargout == 0
    loadMovie = false; loadHis = false; loadProjs = false;
elseif nargout == 1
    loadHis = false; loadProjs = false;
elseif nargout == 2
    loadProjs = false;
end

loadProjs = false; makeProjs = false;

thisExperiment = liveExperiment(Prefix);

[xSize, ySize, ~, ~, ~,...
    nFrames, nSlices, nDigits] = getFrameInfoParams(FrameInfo);

[~,~,~,~, ~,...
    ~, ~, ~,Channel1,Channel2,~,...
    Channel3]...
= readMovieDatabase(Prefix, varargin);
Channels = {Channel1, Channel2, Channel3};

movieMat = []; hisMat = []; maxMat = [];  medMat = []; midMat = [];

nChDatabase = sum(~cellfun(@isempty, Channels)); %this method fails if your
% exported channels don't match your moviedatabase.

nChTifs = 0;
for i = 1:3
    nChTifs = nChTifs + ~isempty(dir([PreProcPath, filesep, Prefix,filesep,'*ch0',num2str(i),'*.tif']));
end

nCh = max(nChTifs, nChDatabase);

pth = [PreProcPath, filesep, Prefix, filesep,Prefix];

tic

if loadMovie && exist([pth, '_movieMat.Mat'], 'file')
    disp('Loading movie mats...')
       movieMat = getMovieMat(thisExperiment);
    disp(['Movie mats loaded. ', num2str(toc), ' s elapsed.'])
    if isempty(movieMat)
        makeMovie = true;
    end
end

if loadHis && exist([pth, '_hisMat.Mat'], 'file')
    
    disp('Loading nuclear mats...')

    
    hisMat = getHisMat(thisExperiment);
    disp(['Nuclear mats loaded. ', num2str(toc), ' s elapsed.'])
    
    if isempty(hisMat)
        makeHis = true;
    end
    
end


if makeMovie
    
    disp('Creating movie mats...')
    
%     startParallelPool(nWorkers, 0, 1);
    
    movieMat = zeros(ySize, xSize,nSlices+nPadding, nFrames, nCh, 'uint16'); % y x z t ch
    
    if makeHis
        hisMat = zeros(ySize, xSize, nFrames, 'uint8'); % y x t
    end
    
    for ch = 1:nCh
        for f = 1:nFrames
            
            for z = 1:nSlices+nPadding
                movieMat(:, :, z, f, ch) = imread([pth,'_',iIndex(f, nDigits),...
                    '_z', iIndex(z, 2), ['_ch', iIndex(ch, 2)], '.tif']);
            end
            
            if makeHis
                hisMat(:, :, f) = imread([pth,'-His_', iIndex(f, nDigits), '.tif']);
            end
            
        end
    end
    
%     save([pth, '_movieMat.Mat'],'movieMat', '-v7.3', '-nocompression');
%     save([pth, '_hisMat.Mat'],'hisMat', '-v7.3', '-nocompression');
    
    disp(['Movie mats created.' , num2str(toc), ' s elapsed.'])
    
end

if  makeHis && ~makeMovie
    
    hisMat = zeros(xSize, ySize,nFrames, 'uint8'); % y x t
    
    for f = 1:nFrames
        hisMat(:, :, f) = imread([pth,'-His_', iIndex(f, nDigits), '.tif']);
    end
    
%     save([pth, '_hisMat.Mat'],'hisMat', '-v7.3', '-nocompression');
    
end

if loadProjs
    
    if exist([pth, '_maxMat.Mat'], 'file') && exist([pth, '_midMat.Mat'], 'file')
        
        load([pth, '_maxMat.Mat'],'maxMat');
        %     load([pth, '_medMat.Mat'],'medMat');
        load([pth, '_midMat.Mat'],'midMat');
        
        if isempty(maxMat) || isempty(midMat)
            makeProjs = true;
        end
        
    end
    
    if makeProjs
%         maxMat = squeeze(max(movieMat(:,:,:,:, :), [], 3)); % y x z t ch
%         midMat = squeeze(max(movieMat(:,:,round(nSlices * .50):round(nSlices * .75),:, :), [], 3));
        save([pth, '_maxMat.Mat'],'maxMat', '-v7.3', '-nocompression');
        save([pth, '_midMat.Mat'],'maxMat', '-v7.3', '-nocompression');
    end
    
end


end