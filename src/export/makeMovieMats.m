function [movieMat, hisMat]...
    = makeMovieMats(Prefix, ~, ~, ~, varargin)

makeHis = false;
makeMovie = false;

loadHis = true;
loadMovie = true;

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
    loadMovie = false; loadHis = false;
elseif nargout == 1
    loadHis = false;
end

thisExperiment = liveExperiment(Prefix);
FrameInfo = getFrameInfo(thisExperiment);
preFolder = thisExperiment.preFolder;

[xSize, ySize, ~, ~, ~,...
    nFrames, nSlices, nDigits] = getFrameInfoParams(FrameInfo);

Channels = thisExperiment.getChannels();

movieMat = []; hisMat = [];

nChDatabase = sum(~cellfun(@isempty, Channels)); %this method fails if your
% exported channels don't match your moviedatabase.

nChTifs = 0;
for i = 1:3
    nChTifs = nChTifs +...
        ~isempty(dir([preFolder, filesep,Prefix,'*ch0',num2str(i),'*.tif']));
end

nCh = max(nChTifs, nChDatabase);

tic

if loadMovie && exist([preFolder,filesep,Prefix, '_movieMat.Mat'], 'file')
    disp('Loading movie mats...')
    movieMat = getMovieMat(thisExperiment);
    disp(['Movie mats loaded. ', num2str(toc), ' s elapsed.'])
    if isempty(movieMat)
        makeMovie = true;
    end
end

if loadHis && exist([preFolder,filesep, Prefix, '_hisMat.Mat'], 'file')
    
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
        parfor f = 1:nFrames
            
            for z = 1:nSlices+nPadding
                movieMat(:, :, z, f, ch) = imread([preFolder, filesep, Prefix, '_',iIndex(f, nDigits),...
                    '_z', iIndex(z, 2), ['_ch', iIndex(ch, 2)], '.tif']);
            end
            
            if makeHis
                hisMat(:, :, f) = imread([preFolder, filesep, Prefix, '-His_', iIndex(f, nDigits), '.tif']);
            end
            
        end
    end
    
    livemRNAImageMatSaver([preFolder, filesep, Prefix, '_movieMatCh1.mat'],...
        movieMat(:, :, :, :, 1));
    if size(movieMat, 5) > 1
        livemRNAImageMatSaver([preFolder, filesep, Prefix, '_movieMatCh2.mat'],...
            movieMat(:, :, :, :, 2));
    end
    if size(movieMat, 5) == 3
        livemRNAImageMatSaver([preFolder, filesep, Prefix, '_movieMatCh3.mat'],...
            movieMat(:, :, :, :, 3));
    end
    
    livemRNAImageMatSaver([preFolder,filesep, Prefix, '_hisMat.mat'], hisMat);
    
    disp(['Movie mats created.' , num2str(toc), ' s elapsed.'])
    
end

if  makeHis && ~makeMovie
    
    hisMat = zeros(ySize, xSize,nFrames, 'uint8'); % y x t
    
    for f = 1:nFrames
        hisMat(:, :, f) = imread([preFolder, filesep, Prefix, '-His_', iIndex(f, nDigits), '.tif']);
    end
    
    livemRNAImageMatSaver([preFolder,filesep, Prefix, '_hisMat.mat'], hisMat);
    
end

end