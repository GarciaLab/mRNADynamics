function [FrameInfo, PixelSize] = loadFrameInfo(OutputFolder, PreProcPath, Prefix)

FrameInfoPath = [OutputFolder, filesep, 'FrameInfo.mat'];

if exist(FrameInfoPath, 'file')
    load(FrameInfoPath, 'FrameInfo');
    
    %See if this came from the 2-photon, which is the default
    if ~isfield(FrameInfo, 'FileMode') || strcmp(FrameInfo(end).FileMode, 'TIF')
        PixelSize = 0.22; %um
    elseif strcmp(FrameInfo(1).FileMode, 'LSM') | strcmp(FrameInfo(1).FileMode, 'LSMExport')%#ok<*OR2>
        PixelSize = FrameInfo(1).PixelSize;
    elseif strcmp(FrameInfo(1).FileMode, 'OMETIFF') || strcmp(FrameInfo(1).FileMode, 'LIFExport') || strcmp(FrameInfo(1).FileMode, 'LAT') || strcmp(FrameInfo(1).FileMode, 'DSPIN')%CS20170907
        PixelSize = FrameInfo(1).PixelSize;
    else
        error('FileMode %s not supported', FrameInfo(1).FileMode);
    end
    
else
    warning('No FrameInfo.mat detected. Trying to pull out magnification information from the TIF file')
    
    DZoom = dir([PreProcPath, filesep, Prefix, filesep, '*z*.tif']);
    ImageInfo = imfinfo([PreProcPath, filesep, Prefix, filesep, DZoom(1).name]);
    PixelSize = 1 / ImageInfo.XResolution;
end

end