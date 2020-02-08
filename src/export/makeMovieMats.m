function [movieMat, hisMat, maxMat, medMat, midMat]...
    = makeMovieMats(Prefix, PreProcPath, nWorkers, FrameInfo, Channels, varargin)

noLoad = false;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'noLoad')
        noLoad = true;
    end
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

nPadding = 2; %normally we pad a blank image above and below the stack.
%if a movie is unpadded this will require modification

pth = [PreProcPath, filesep, Prefix, filesep,Prefix];

movie = false;
maxmat = false;
his = false;

if exist([pth, '_movieMat.Mat'], 'file')
    disp('Loading movie mats...')
    tic
    if ~noLoad
        load([pth, '_movieMat.Mat'],'movieMat');
    end
    disp('Movie mats loaded.')
    toc
    if ~isempty(movieMat)
        movie = true;
    end
end

if exist([pth, '_hisMat.Mat'], 'file')
    if ~noLoad
        load([pth, '_hisMat.Mat'],'hisMat');
    end
    if ~isempty(hisMat)
        his = true;
    end
end


if ~movie
    
    disp('Creating movie mats...')
    tic
    
    startParallelPool(nWorkers, 0, 1);

    movieMat = zeros(nCh, nSlices+nPadding, nFrames, xSize, ySize, 'uint16'); % ch z t x y
    maxMat = zeros(nCh, nFrames, xSize, ySize, 'uint16'); % ch z x y
    maxMat = zeros(nCh, nFrames, xSize, ySize, 'uint16'); % ch z x y
    medMat = zeros(nCh, nFrames, xSize, ySize, 'uint16'); % ch z x y
    hisMat = zeros(nFrames, xSize, ySize, 'uint16'); % f x y
    for ch = 1:nCh
        for f = 1:nFrames
            
            for z = 1:nSlices+nPadding
                movieMat(ch, z, f, :, :) = imread([pth,'_',iIndex(f, nDigits), '_z', iIndex(z, 2), ['_ch', iIndex(ch, 2)], '.tif']);
            end
            
            hisMat(f, :, :) = imread([pth,'-His_', iIndex(f, nDigits), '.tif']);
            
        end
    end
    
    save([pth, '_movieMat.Mat'],'movieMat', '-v7.3', '-nocompression');
    save([pth, '_hisMat.Mat'],'hisMat', '-v7.3', '-nocompression');
    
    disp('Movie mats created.')
    toc
    
end

if movie && ~his
    hisMat = zeros(nFrames, xSize, ySize, 'uint16'); % f x y
    
    parfor f = 1:nFrames
        hisMat(f, :, :) = imread([pth,'-His_', iIndex(f, nDigits), '.tif']);
    end
    
    save([pth, '_hisMat.Mat'],'hisMat', '-v7.3', '-nocompression');
    
end

if exist([pth, '_maxMat.Mat'], 'file') && exist([pth, '_medMat.Mat'], 'file') && exist([pth, '_midMat.Mat'], 'file')
    
    load([pth, '_maxMat.Mat'],'maxMat');
    load([pth, '_medMat.Mat'],'medMat');
    load([pth, '_midMat.Mat'],'midMat');
    
    if ~isempty(maxMat)
        maxmat = true;
    end
    
end

if ~maxmat
    maxMat = squeeze(max(movieMat(:,:,:,:, :), [], 2)); % ch z t x y
    medMat = []; %median is so slow. %     medMat = squeeze(median(movieMat(:,:,:,:, :), 2));
    midMat = squeeze(max(movieMat(:,round(nSlices * .50):round(nSlices * .75),:,:, :), [], 2));
    save([pth, '_maxMat.Mat'],'maxMat', '-v7.3', '-nocompression');
    save([pth, '_medMat.Mat'],'maxMat', '-v7.3', '-nocompression');
    save([pth, '_midMat.Mat'],'maxMat', '-v7.3', '-nocompression');
end

toc

end