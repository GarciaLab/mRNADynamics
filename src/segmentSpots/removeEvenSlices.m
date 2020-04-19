function removeEvenSlices(Prefix, varargin)


[~, ProcPath] = DetermineLocalFolders(Prefix);

% Get all files in the current folder
dogFolder = [ProcPath, filesep, Prefix, '_', filesep, 'dogs', filesep];
files = dir([dogFolder, filesep, '*.tif';]);

for id = 1:length(files)
    
    imFile = files(id).name;
    
    dogStackOld = imreadStack([dogFolder, filesep, imFile]);
    
    ySize = size(dogStackOld, 1);
    xSize = size(dogStackOld, 2);
    zSize = size(dogStackOld, 3);
    
    
    dogStack = zeros(ySize, xSize, zSize/2, 'like', dogStackOld);
    n = 0;
    for zOdd = 1:2:(zSize - 1)
        
        n = n + 1;
        dogStack(:, :, n) = dogStackOld(:, :, zOdd);
        dogStack = uint16(dogStack*10000);
        
        if n == 1
            imwrite(dogStack, [dogFolder,filesep, imFile], 'WriteMode', 'overwrite');
            
        else
            imwrite(dogStack, [dogFolder,filesep, imFile], 'WriteMode','append');
        end
        
        
    end
    
    
end