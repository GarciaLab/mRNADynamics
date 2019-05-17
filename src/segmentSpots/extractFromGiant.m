function dogs = extractFromGiant(giantIm, format, padSize, firstFrame, lastFrame, Prefix, channel, outPath, noSave, varargin)

dim = length(size(giantIm));
frameInterval = padSize+format(2);
dogs = zeros(format(1), format(2), format(3)-2, lastFrame - firstFrame + 1);
padZ = false;

if ~isempty(varargin)
    probs = varargin{1};
end
fcnt = 1;
for frame = firstFrame:lastFrame
    ind1 = frameInterval*(fcnt-1) + 1;
    ind2 = frameInterval + (frameInterval*(fcnt-1)) - padSize;
    
    if dim == 2
        im = giantIm(:,ind1:ind2);
    elseif dim == 3
        
        if probs
            im = gather(uint16(giantIm(:,ind1:ind2, :)));
%             imshow(im(:,:,5),[median(im(:)),max(im(:))]);
        else
            im = gather(uint16((giantIm(:,ind1:ind2, :) + 100)*10));
        end
        
        %pad the image to avoid edge effects
        im(:, 1:1+round((padSize/2)), :) = 0;
        im(:, end-round((padSize/2)):end, :) = 0;
        %         imshow(im(:, :, 6), [median(im(:)), max(im(:))]);
        
        if padZ
            im = cat(zeros(size(im,1), size(im,2)), im, 3);
            im(:,:,end+1) = zeros(size(im,1), size(im,2));
        end
        
        if noSave
            dogs(:, :, :, fcnt) = im;
        end
        if probs
            fldr = 'custProbs';
        else
            fldr = 'dogs';
        end
        nameSuffix = ['_ch', iIndex(channel, 2)];
        
        mkdir([outPath, filesep,Prefix,'_',filesep,fldr]);
        
        
        for z = 1:size(im,3)
            plane = im(:,:,z);
            if ~noSave
                
                dog_name = ['DOG_', Prefix, '_', iIndex(frame, 3), '_z', iIndex(z, 2), nameSuffix, '.tif'];
                dog_full_path = [outPath, filesep,Prefix,'_',filesep,fldr,filesep,dog_name];
                %         imwrite(uint16(dog(:,:, z)), dog_full_path);
                imwrite(plane,dog_full_path);
                % save(['E:\Armando\giantDogs\','dog_frame',num2str(frame),'_z',num2str(z),'.mat'], 'plane');
            end
        end
    end
    
    fcnt = fcnt + 1;
    
end

end