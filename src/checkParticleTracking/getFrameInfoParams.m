function [xSize_px, ySize_px, pixelSize_nm, zStep_um, snippet_size_px,...
    nFrames, nSlices, nDigits] = getFrameInfoParams(FrameInfo)
  
  xSize_px = mean([FrameInfo.PixelsPerLine]);
  ySize_px = mean([FrameInfo(1).LinesPerFrame]);
  pixelSize_nm = mean([FrameInfo.PixelSize]) * 1000; %nm
  zStep_um = mean([FrameInfo(1).ZStep]);
  nSlices = mean([FrameInfo(1).NumberSlices]);
  snippet_size_px = 2 * (floor(1300 / (2 * pixelSize_nm))) + 1; % nm. note that this is forced to be odd
  nFrames = length(FrameInfo);
  
  
%See how  many frames we have and adjust the index size of the files to load accordingly
if nFrames < 1E3
    nDigits = 3;
elseif nFrames < 1E4
    nDigits = 4;
else
    error('No more than 10,000 frames supported.')
end
  
end