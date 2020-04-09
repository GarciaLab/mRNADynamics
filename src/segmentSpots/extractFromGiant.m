function dwarfImage = extractFromGiant(giantImage, format, padSize, frameRange, ...
    Prefix, channel, outPath, noSave, varargin)

dim = length(size(giantImage));
frameInterval = padSize+format(2);
firstFrame = frameRange(1);
lastFrame = frameRange(2);
dwarfImage = zeros(format(1), format(2), format(3)-2, lastFrame - firstFrame + 1);
padZ = false;
probs = false;
numType = 'single';
mat = false;
% saveType = '.tif';
saveType = '5D';

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'single')
        numType = 'single';
    elseif strcmpi(varargin{i}, 'double')
        numType = 'double';
    elseif strcmpi(varargin{i}, 'probs')
        probs = true;
    elseif strcmpi(varargin{i}, 'mat')
        mat = true;
        saveType = '.mat';
    elseif strcmpi(varargin{i}, 'padZ')
        padZ = true;
    elseif strcmpi(varargin{i}, 'noSave')
        noSave = true;
        saveType = 'none';
    end
end


fcnt = 1;
for frame = firstFrame:lastFrame
    ind1 = frameInterval*(fcnt-1) + 1;
    ind2 = frameInterval + (frameInterval*(fcnt-1)) - padSize;
    
    if dim == 2
        im = giantImage(:,ind1:ind2);
    elseif dim == 3
        
        switch saveType
            case '.tif'
                if probs
                    im = gather(uint16(giantImage(:,ind1:ind2, :)));
                    %             imshow(im(:,:,5),[median(im(:)),max(im(:))]);
                else
                    im = gather(uint16((giantImage(:,ind1:ind2, :) + 100)*10));
                end
            case {'.mat', 'none'}
                im = gather(giantImage(:,ind1:ind2, :));
        end
        
        
        
        %pad the image to avoid edge effects
        im(:, 1:1+round((padSize/2)), :) = 0;
        im(:, end-round((padSize/2)):end, :) = 0;
        %         imshow(im(:, :, 6), [median(im(:)), max(im(:))]);
        
        if padZ
            im = cat(zeros(size(im,1), size(im,2)), im, 3);
            im(:,:,format(3)) = zeros(size(im,1), size(im,2));
        end
        
        dwarfImage(:, :, :, fcnt) = im;
        
        if probs
            fldr = 'custProbs';
        else
            fldr = 'dogs';
        end
        nameSuffix = ['_ch', iIndex(channel, 2)];
        
        mkdir([outPath, filesep,Prefix,'_',filesep,fldr]);
        
        if strcmpi(saveType, '5D')
            save([outPath, filesep,Prefix,'_',filesep,'_dogMat.mat'], dwarfImage, '-v7.3', '-nocompression');
        else
            
            for z = 1:size(im,3)
                
                plane = im(:,:,z);
                dog_name = ['DOG_', Prefix, '_', iIndex(frame, 3), '_z', iIndex(z, 2), nameSuffix];
                dog_full_path = [outPath, filesep,Prefix,'_',filesep,fldr,filesep,dog_name];
                
                switch saveType
                    case '.tif'
                        imwrite(plane,[dog_full_path, saveType]);
                    case '.mat'
                        save([dog_full_path,saveType], 'plane');
                    case 'none'
                        %do nothing
                end
                
            end
        end
        
    end
    
    
    fcnt = fcnt + 1;
    
end

end
