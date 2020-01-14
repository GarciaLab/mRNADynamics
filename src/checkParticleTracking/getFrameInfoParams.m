function [xSize, ySize, pixelSize, zStep, snippet_size, LinesPerFrame, PixelsPerLine,...
    numFrames] = getFrameInfoParams(FrameInfo)
  
  xSize = FrameInfo(1).PixelsPerLine;
  ySize = FrameInfo(1).LinesPerFrame;
  pixelSize = FrameInfo(1).PixelSize * 1000; %nm
  zStep = FrameInfo(1).ZStep;
  snippet_size = 2 * (floor(1300 / (2 * pixelSize))) + 1; % nm. note that this is forced to be odd
  LinesPerFrame = FrameInfo(1).LinesPerFrame;
  PixelsPerLine = FrameInfo(1).PixelsPerLine;
  numFrames = length(FrameInfo);

end