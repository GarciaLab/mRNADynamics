function dogMat = loadDogMat(inputString, varargin)

warning('off', 'MATLAB:MatFile:OlderFormat')

disp(['Loading dog movie: ',inputString,'...']);
tic;

frameRange = [];
isWritable = false;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end


%either pass the moviefile path directly or
% load with the project prefix
if ~contains(inputString, '.mat')
    Prefix = inputString;
    [~, ProcPath] = DetermineLocalFolders(Prefix);
    dogFile = [ProcPath, filesep, Prefix, '_', filesep, Prefix, '_dogMat.mat'];
else
    dogFile = inputString;
end

movieMat = loadMovieMat(Prefix);
ySize = size(movieMat, 1);
xSize = size(movieMat, 2);
nFrames = size(movieMat, 4);
nSlices =  size(movieMat, 3);
clear movieMat;
ch = 2;
if ~exist(dogFile, 'file')
    dogMat = [];
    for f = 1:nFrames
        for z = 1:nSlices
            load(...
                [ProcPath, filesep, filesep, Prefix, '_', filesep, 'dogs', filesep,...
                'DOG_', Prefix,'_',iIndex(f, 3), '_z', iIndex(z, 2)...
                '_ch', iIndex(ch, 2),'.mat'], 'plane');
            dogMat(:, :, z, f) = plane;
        end
    end
    
    if whos(var2str(dogMat)).bytes < 2E9
        save(dogFile, 'dogMat', '-v6');
    else
        dogMatic = newmatic(dogFile,true,...
            newmatic_variable('dogMat', 'double', [ySize, xSize, nSlices, nFrames], [ySize, xSize, 1, 1]));
        dogMatic.dogMat = dogMat;
    end
    return
    
else
    
    dogmatfile = matfile(dogFile, 'Writable', isWritable);
    dims = size(dogmatfile, 'dogMat');
    
end


if isempty(frameRange)
    frameRange = 1:dims(3);
else
    frameRange = frameRange(1):frameRange(end);
end

dogMat = dogmatfile.dogMat(:, :, frameRange);

disp(['Nuclear movie loaded- ', num2str(toc), 's']);


end
