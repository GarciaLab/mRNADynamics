function hisMat = loadHisMat(hisFile,  varargin)

disp('Loading nuclear movie....');

frameRange = [];
isWritable = false;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

hismatfile = matfile(hisFile, 'Writable', isWritable);
dims = size(hismatfile.hisMat);


if isempty(frameRange)
    frameRange = 1:dims(3);
else
    frameRange = frameRange(1):frameRange(end);
end

hisMat = hismatfile.hisMat(:, :, frameRange);

disp('Nuclear movie loaded.');


end