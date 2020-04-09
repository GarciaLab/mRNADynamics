function format = preProcFormat(PreProcFolder, Prefix)

format = {};
pth = [PreProcFolder, filesep, Prefix];

%possible file types

mat2D = ~isempty(dir([pth, filesep,'*_z*.mat'])); if mat2D, format = [format var2str(mat2D)]; end
mat3D = ~isempty([]); if mat3D, format = [format var2str(mat3D)]; end %reserving this for the future
mat5D = ~isempty(dir([pth, filesep,'*movieMat*.mat']));  if mat5D, format = [format var2str(mat5D)]; end

tif2D = ~isempty(dir([pth, filesep,'*_z*.tif'])); if tif2D, format = [format var2str(tif2D)]; end
tif3D = ~isempty(dir([pth, filesep, 'stacks', filesep, '*.tif'])); if tif3D, format = [format var2str(tif3D)]; end

if length(format) == 1
    format = format{1};
end

end