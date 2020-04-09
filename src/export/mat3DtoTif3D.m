function mat3DtoTif3D(mat3D, PreProcFolder, Prefix)

for f = 1:size(mat3D, 1)
    imwrite(uint16(squeeze(mat3D(f, :, :))), [PreProcFolder, filesep, Prefix, filesep, Prefix, '_hisMovie.tif'], 'WriteMode', 'append')
end

end