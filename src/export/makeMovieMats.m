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

[xSize, ySize, ~, ~, ~,...
    nFrames, nSlices, nDigits] = getFrameInfoParams(FrameInfo);


movieMat = []; hisMat = []; maxMat = [];  medMat = []; midMat = [];

% nCh = sum(~cellfun(@isempty, Channels)); %this method fails if your
% exported channels don't match your moviedatabase.

nCh = 0;
for i = 1:3
    nCh = nCh + ~isempty(dir([PreProcPath, filesep, Prefix,filesep,'*ch0',num2str(i),'*.tif']));
end

pth = [PreProcPath, filesep, Prefix, filesep,Prefix];

tic

if loadMovie && exist([pth, '_movieMat.Mat'], 'file')
    disp('Loading movie mats...')
        load([pth, '_movieMat.Mat'],'movieMat');
    disp(['Movie mats loaded. ', num2str(toc), ' s elapsed.'])
    if isempty(movieMat)
        makeMovie = true;
    end
end

if loadHis && exist([pth, '_hisMat.Mat'], 'file')
    
    disp('Loading nuclear mats...')

    load([pth, '_hisMat.Mat'],'hisMat');
    disp(['Nuclear mats loaded. ', num2str(toc), ' s elapsed.'])
    
    if isempty(hisMat)
        makeHis = true;
    end
    
end


if makeMovie
    
    disp('Creating movie mats...')
    
    startParallelPool(nWorkers, 0, 1);
    
    movieMat = zeros(nCh, nSlices+nPadding, nFrames, ySize, xSize, 'uint16'); % ch z t x y
    
    if makeHis
        hisMat = zeros(nFrames, ySize, xSize, 'uint16'); % f x y
    end
    
    for ch = 1:nCh
        for f = 1:nFrames
            
            parfor z = 1:nSlices+nPadding
                movieMat(ch, z, f, :, :) = imread([pth,'_',iIndex(f, nDigits), '_z', iIndex(z, 2), ['_ch', iIndex(ch, 2)], '.tif']);
            end
            
            if makeHis
                hisMat(f, :, :) = imread([pth,'-His_', iIndex(f, nDigits), '.tif']);
            end
            
        end
    end
    
    save([pth, '_movieMat.Mat'],'movieMat', '-v7.3', '-nocompression');
    save([pth, '_hisMat.Mat'],'hisMat', '-v7.3', '-nocompression');
    
    disp(['Movie mats created.' , num2str(toc), ' s elapsed.'])
    
end

if  makeHis && ~makeMovie
    
    hisMat = zeros(nFrames, xSize, ySize, 'uint16'); % f x y
    
    parfor f = 1:nFrames
        hisMat(f, :, :) = imread([pth,'-His_', iIndex(f, nDigits), '.tif']);
    end
    
    save([pth, '_hisMat.Mat'],'hisMat', '-v7.3', '-nocompression');
    
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
        maxMat = squeeze(max(movieMat(:,:,:,:, :), [], 2)); % ch z t x y
        %     medMat = []; %median is so slow. %     medMat = squeeze(median(movieMat(:,:,:,:, :), 2));
        midMat = squeeze(max(movieMat(:,round(nSlices * .50):round(nSlices * .75),:,:, :), [], 2));
        save([pth, '_maxMat.Mat'],'maxMat', '-v7.3', '-nocompression');
        %     save([pth, '_medMat.Mat'],'maxMat', '-v7.3', '-nocompression');
        save([pth, '_midMat.Mat'],'maxMat', '-v7.3', '-nocompression');
    end
    
end


end