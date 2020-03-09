function hisMat = loadHisMat(hisFile,  dims, varargin)

frameRange = [];

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

ySize = dims(1);
xSize = dims(2);
nFrames = dims(4);

hisMatic = newmatic(hisFile,...
            newmatic_variable('hisMat', 'uint16', [ySize, xSize, nFrames], [ySize, xSize, 1]));

    if ~isempty(frameRange)
        hisMat = hisMatic.hisMat(:, :, frameRange);
    else
        hisMat = hisMatic.hisMat;
    end

end