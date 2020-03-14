function movieMat = loadMovieMat(inputString, varargin)

warning('off', 'MATLAB:MatFile:OlderFormat')

disp(['Loading movie: ',inputString,'...']);
tic;

zRange = [];
frameRange = [];
chRange = [];
movieMatCh1=[]; movieMatCh2=[]; movieMatCh3=[];
isWritable = false;
isDividedIntoChannels = false;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end


%either pass the moviefile path directly or...
% load with the project prefix
if ~contains(inputString, '.mat')
    
    Prefix = inputString;
    [~, ~, ~, ~, PreProcPath] = DetermineLocalFolders(Prefix);
    
    PreProcFolder = [PreProcPath, filesep, Prefix, filesep];
    
    movieChDir = dir([PreProcFolder, filesep, Prefix, '_movieMatCh*.mat']);
    numChannelsToLoad = numel(movieChDir);
    
    if ~isempty(chRange) && max(chRange) < 3
        numChannelsToLoad = numChannelsToLoad-1;
    end
    
    if numChannelsToLoad > 0
        isDividedIntoChannels = true;

        for ch = 1:numChannelsToLoad
            load([PreProcFolder, filesep, movieChDir(ch).name]);
        end
        
        try movieMat = squeeze(cat(5, movieMatCh1, movieMatCh2, movieMatCh3));
        catch movieMat =  squeeze(cat(5, movieMatCh1, movieMatCh2)); end
        
        dims = size(movieMat);
        
    else
        movieFile = [PreProcPath, filesep, Prefix, filesep, Prefix, '_movieMat.mat'];
    end
    
else
    movieFile = inputString;
end

if ~isDividedIntoChannels
    moviematfile = matfile(movieFile, 'Writable', isWritable);
    dims = size(moviematfile, 'movieMat');
end


if isempty(frameRange)
    frameRange = 1:dims(4);
else
    frameRange = frameRange(1):frameRange(end);
end

if isempty(zRange)
    zRange = 1:dims(3);
else
    zRange = zRange(1):zRange(end);
end

if isempty(chRange)
    chRange = 1:dims(5);
else
    chRange = chRange(1):chRange(end);
end

if isDividedIntoChannels
    movieMat = movieMat(:, :, zRange, frameRange, chRange);
else
    movieMat = moviematfile.movieMat(:, :, zRange, frameRange, chRange);
end

disp(['Movie loaded- ', num2str(toc), 's']);


end
