function [xSize, ySize, pixelSize, zStep, snippet_size, LinesPerFrame, PixelsPerLine,...
    numFrames, NDigits, NChannels, Particles, Spots, SpotFilter] = getFrameInfoParams(FrameInfo, Particles, Spots, SpotFilter)
  
  xSize = FrameInfo(1).PixelsPerLine;
  ySize = FrameInfo(1).LinesPerFrame;
  pixelSize = FrameInfo(1).PixelSize * 1000; %nm
  zStep = FrameInfo(1).ZStep;
  snippet_size = 2 * (floor(1300 / (2 * pixelSize))) + 1; % nm. note that this is forced to be odd
  LinesPerFrame = FrameInfo(1).LinesPerFrame;
  PixelsPerLine = FrameInfo(1).PixelsPerLine;
  numFrames = length(FrameInfo);
  
  %See how  many frames we have and adjust the index size of the files to load accordingly
  if numFrames < 1E3
    NDigits = 3;
  elseif numFrames < 1E4
    NDigits = 4;
  else
    error('No more than 10,000 frames supported.')
  end
  
  %Create the particle array. This is done so that we can support multiple
  %channels. Also figure out the number of channels
  if iscell(Particles)
    NChannels = length(Particles);
  else
    Particles = {Particles};

    if ~iscell(Spots)
      Spots = {Spots};
    end

    SpotFilter = {SpotFilter};
    NChannels = 1;
  end
end