function removeEvenSlices(Prefix, varargin)

saveAsMat = false;
saveAsStack = false;

for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

[~, ProcPath] = DetermineLocalFolders(Prefix);

% Get all files in the current folder
dogFolder = [ProcPath, filesep, Prefix, '_', filesep, 'dogs', filesep];
files = dir([dogFolder, filesep, '*.tif';]);

for id = 1:length(files)
    
    imFile = files(id).name;
    
    dogStackOld = readTiffStack([dogFolder, filesep, imFile]);
    
    ySize = size(dogStackOld, 1);
    xSize = size(dogStackOld, 2);
    zSize = size(dogStackOld, 3);
    
    dogStack = zeros(ySize, xSize, zSize/2, 'like', dogStackOld);
    n = 0;
    for zOdd = 1:2:(zSize - 1)
        n = n + 1;
        dogStack(:, :, n) = dogStackOld(:, :, zOdd);
        
        if ~saveAsStack && ~saveAsMat
            dog = dogStack(:, :, n);
            channelNameCorrection = strrep(imFile, '_ch', ['_z', iIndex(n, 2),'_ch']);
            imName = ['prob', Prefix, '_' channelNameCorrection];
			imwrite(uint16(dog*10000), [dogFolder,filesep, imName, '.tif']);
        end
    end
    
    if saveAsMat
        save([dogFolder,filesep, imFile, '.mat'], 'dogStack', '-v6');
    end    
end