function removeEvenSlices(Prefix, varargin)

saveType = 'mat3D';

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
    
    dogStackOld = imreadStack([dogFolder, filesep, imFile]);
    
    ySize = size(dogStackOld, 1);
    xSize = size(dogStackOld, 2);
    zSize = size(dogStackOld, 3);
    
    
    dogStack = zeros(ySize, xSize, zSize/2, 'like', dogStackOld);
    n = 0;
    for zOdd = 1:2:(zSize - 1)
        
        n = n + 1;
        dogStack(:, :, n) = dogStackOld(:, :, zOdd);
        
        %if we're saving as individual planes, do this inside the loop
        if ~strcmpi(saveType, 'mat3D')
            
            imName = strrep(imFile(1:end-4), '_ch', ['_z', iIndex(n, 2),'_ch']);
            if n == 1
                imwrite(uint16(dogStack*10000), [dogFolder,filesep, imName, '.tif']);
            else
                imwrite(uint16(dogStack*10000), [dogFolder,filesep, imName, '.tif'], 'WriteMode','append');
                delete([dogFolder, filesep, imFile]);
            end
            
        end
        
    end %loop over planes
    
    
    %if we're saving as a stack, save here
    if strcmpi(saveType, 'mat3D')
        dogStack = uint16(dogStack*10000);
        save([dogFolder,filesep, imFile(1:end-4), '.mat'], 'dogStack', '-v6');
    end
    
    
end %loop over frames



end