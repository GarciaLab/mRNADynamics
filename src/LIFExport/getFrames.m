function [NSeries, NFrames, NSlices, NPlanes, NChannels, Frame_Times] = getFrames(LIFMeta)
  
  NSeries = LIFMeta.getImageCount();
  
  %Figure out the number of slices in each series
  NSlices = [];
  for i = 1:NSeries
      NSlices(i) = str2num(LIFMeta.getPixelsSizeZ(i-1));
  end

  %Number of planes per series
  NPlanes = [];
  for i = 1:NSeries
      NPlanes(i) = LIFMeta.getPlaneCount(i-1);
  end

  %Number of channels
  NChannels = LIFMeta.getChannelCount(0);

  %Finally, use this information to determine the number of frames in each series        
  NFrames = NPlanes./NSlices/NChannels;

  %Get rid of the last frame as it is always incomplete because that's when we stopped it
  NFrames = NFrames - 1;
  NPlanes = NPlanes - NSlices * NChannels;      
  Frame_Times = zeros(1, sum(NFrames.*NSlices));
end
