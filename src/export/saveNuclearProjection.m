function saveNuclearProjection(hisMat, file)

nFrames = size(hisMat, 3);

for f = 1:nFrames
    
    image = hisMat(:, :, f);
    if f == 1
        imwrite(image,file, 'Compression', 'none');
    else
        imwrite(image,file,'WriteMode', 'append', 'Compression', 'none');
    end
     
end

end