function movieMat = loadMovieMat(inputString, varargin)

warning('off', 'MATLAB:MatFile:OlderFormat')

disp(['Loading movie: ',inputString,'...']);
tic;

zRange = [];
frameRange = [];
chRange = [];
isWritable = false;

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
    movieFile = [PreProcPath, filesep, Prefix, filesep, Prefix, '_movieMat.mat'];
else
    movieFile = inputString;
end
    

moviematfile = matfile(movieFile, 'Writable', isWritable);
dims = size(moviematfile, 'movieMat');


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

movieMat = moviematfile.movieMat(:, :, zRange, frameRange, chRange);

disp(['Movie loaded- ', num2str(toc), 's']);


end
