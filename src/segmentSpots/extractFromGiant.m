function dogs = extractFromGiant(giantIm, format, padSize, firstFrame, lastFrame, Prefix, channel, outPath, noSave)

dim = length(size(giantIm));
frameInterval = padSize+format(2);
dogs = zeros(format(1), format(2), format(3)-2, lastFrame - firstFrame + 1);

fcnt = 1;
for frame = firstFrame:lastFrame
    ind1 = frameInterval*(fcnt-1) + 1;
    ind2 = frameInterval + (frameInterval*(fcnt-1)) - padSize;
    
    if dim == 2
        im = giantIm(:,ind1:ind2);
    elseif dim == 3
        
        im = gather(uint16((giantIm(:,ind1:ind2, :) + 100)*10));
        im(:, end - (padSize/2) : end, :) = 0;
        if noSave
            dogs(:, :, :, fcnt) = im; 
        end
        
        for z = 1:size(im,3)
            plane = im(:,:,z);
%             imshow(plane,[median(plane(:)), max(plane(:))]);
            if ~noSave
                nameSuffix = ['_ch', iIndex(channel, 2)];
                dog_name = ['DOG_', Prefix, '_', iIndex(frame, 3), '_z', iIndex(z, 2), nameSuffix, '.tif'];
                dog_full_path = [outPath, filesep,Prefix,'_',filesep,'dogs',filesep,dog_name];
                %         imwrite(uint16(dog(:,:, z)), dog_full_path);
                imwrite(plane,dog_full_path)
                % save(['E:\Armando\giantDogs\','dog_frame',num2str(frame),'_z',num2str(z),'.mat'], 'plane');
            end
        end
    end
    
    fcnt = fcnt + 1;
    
end
end