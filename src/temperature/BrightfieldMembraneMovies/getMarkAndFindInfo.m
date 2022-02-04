function [NSeries, NSlices, NPlanes, NChannels, NReplicates] = getMarkAndFindInfo(LIFMeta, ExperimentType)
    if ~exist('ExperimentType', 'var')
        ExperimentType =  {''};
    end
  
  NSeries = LIFMeta.getImageCount();
  
  
  %Figure out the number of slices in each series
  NSlices = [];
  for i = 1:NSeries
      NSlices(i) = str2double(LIFMeta.getPixelsSizeZ(i-1));
  end
  
    NChannels = LIFMeta.getChannelCount(0);


  %Number of planes per series
  NPlanes = [];
  for i = 1:NSeries
      NPlanes(i) = LIFMeta.getPlaneCount(i-1);
  end

  NReplicates = NPlanes./NSlices/NChannels;

end
