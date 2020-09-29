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

  %MT 2020-09-29: modification to process Modular Enhancer 45uW data, which
  %was taken with a fixed number of time frames (2)
  if NSeries == 1 && NFrames == 2
       %NFrames stays the same
      disp('Keeping the last frame of the only series in this dataset.')
  else
    %Get rid of the last frame as it is always incomplete because we stop
    %it in middle of the last frame
    NFrames = NFrames - 1;
    NPlanes = NPlanes - (NSlices * NChannels);
  end
  
  Frame_Times = zeros(1, sum(NFrames.*NSlices));
end
