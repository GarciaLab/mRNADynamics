function movieMat = loadMovieMat(movieFile, dims, varargin)

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

ySize = dims(1);
xSize = dims(2);
zSize = dims(3);
nFrames = dims(4);
nCh = dims(5);

movieMatic = newmatic(movieFile,...
            newmatic_variable('movieMat', 'uint16', [ySize, xSize, zSize, nFrames,  nCh], [ySize, xSize, 1, 1, 1]));

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

movieMat = movieMatic.movieMat(:, :, zRange, frameRange, chRange);

end
