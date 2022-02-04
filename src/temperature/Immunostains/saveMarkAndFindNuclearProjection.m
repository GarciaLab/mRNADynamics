function saveMarkAndFindNuclearProjection(hisMarkAndFindMat, file)

nFrames = size(hisMarkAndFindMat, 4);
nReplicates = size(hisMarkAndFindMat,3);

if max(hisMarkAndFindMat(:)) < 255
    hisMarkAndFindMat = uint8(hisMarkAndFindMat);
else
    hisMarkAndFindMat = uint16(hisMarkAndFindMat);
end
for r=1:nReplicates
    for f = 1:nFrames
        
        
        image = hisMarkAndFindMat(:, :,r, f);
        if f == 1
            imwrite(image,file);
        else
            imwrite(image,file,'WriteMode', 'append');
            
        end
    end
    
end

end