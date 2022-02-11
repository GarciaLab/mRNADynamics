function [NSeries, NFrames, NSlices, NPlanes, NChannels, Frame_Times] = getFrames(LIFMeta)
  
  NSeries = LIFMeta.getImageCount();
  
  %Figure out the number of slices in each series
  NSlices = [];
  for i = 1:NSeries
      NSlices(i) = str2double(LIFMeta.getPixelsSizeZ(i-1));
  end

  %Number of planes per series
  NPlanes = [];
  for i = 1:NSeries
      NPlanes(i) = LIFMeta.getPlaneCount(i-1);
  end

  NChannels = LIFMeta.getChannelCount(0);

  %Finally, use this information to determine the number of frames in each series        
  NFrames = NPlanes./NSlices/NChannels;

  %Get rid of the last frame as it is always incomplete because that's when
  %we stopped it.
  %AR 2/7/21- sometimes we actually want to keep that frame. Let's assume
  %that if there's only one frame per series, we want to keep it. 
  if max(NFrames) ~= 1
      NFrames = NFrames - 1;
      NPlanes = NPlanes - (NSlices * NChannels);
  end
  Frame_Times = zeros(1, sum(NFrames.*NSlices));
  
end
