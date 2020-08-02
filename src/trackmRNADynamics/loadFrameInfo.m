function [FrameInfo, PixelSize_um] = loadFrameInfo(OutputFolder, PreProcPath, Prefix)

FrameInfoPath = [OutputFolder, filesep, 'FrameInfo.mat'];

if exist(FrameInfoPath, 'file')
    load(FrameInfoPath, 'FrameInfo');
    
    %See if this came from the 2-photon, which is the default
    if ~isfield(FrameInfo, 'FileMode') || strcmp(FrameInfo(end).FileMode, 'TIF')
        PixelSize_um = 0.22; %um
    elseif strcmp(FrameInfo(1).FileMode, 'LSM') | strcmp(FrameInfo(1).FileMode, 'LSMExport')%#ok<*OR2>
        PixelSize_um = FrameInfo(1).PixelSize;
    elseif strcmp(FrameInfo(1).FileMode, 'OMETIFF') || strcmp(FrameInfo(1).FileMode, 'LIFExport') || strcmp(FrameInfo(1).FileMode, 'LAT') || strcmp(FrameInfo(1).FileMode, 'DSPIN')%CS20170907
        PixelSize_um = FrameInfo(1).PixelSize;
    else
        error('FileMode %s not supported', FrameInfo(1).FileMode);
    end
    
else
    warning('No FrameInfo.mat detected. Trying to pull out magnification information from the TIF file')
    
    DZoom = dir([PreProcPath, filesep, Prefix, filesep, '*z*.tif']);
    ImageInfo = imfinfo([PreProcPath, filesep, Prefix, filesep, DZoom(1).name]);
    xRes = [ImageInfo.XResolution];
    PixelSize_um = 1 / xRes(1);
    FrameInfo = [];
end

end