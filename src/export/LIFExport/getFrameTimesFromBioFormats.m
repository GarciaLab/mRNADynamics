function Frame_Times = getFrameTimesFromBioFormats(LIFMeta, NSlices)

import loci.common.DateTools;

nChannels = LIFMeta.getChannelCount(0);

nSeries = LIFMeta.getImageCount();

Frame_Times = [];
movieT0 = nan;

for seriesIndex = 0:nSeries-1
    
    nPlanes = LIFMeta.getPlaneCount(seriesIndex);
    
    %MT 2020-09-29: modification to process Modular Enhancer 45uW data, 
    %which was taken with a fixed number of time frames (2) in a single
    %series
    if nSeries ~= 1
        nPlanes = nPlanes - (NSlices(seriesIndex+1)*nChannels);
    end
    
    if nPlanes > 0
        seriesStamp = LIFMeta.getImageAcquisitionDate(seriesIndex);
        seriesStamp_char = seriesStamp.getValue();
        seriesT0 = DateTools.getTime(seriesStamp_char, DateTools.ISO8601_FORMAT) / 1000; %unix epoch ms -> s
        
        %most of the time this will occur at series 0 but not always
        if isnan(movieT0)
            movieT0 = seriesT0; 
        end
        
        
        for plane = 0:nPlanes-1
            
            deltaT = LIFMeta.getPlaneDeltaT(seriesIndex, plane);
            if ~isempty(deltaT)
                Frame_Times = [Frame_Times, seriesT0 + double(deltaT.value)];
            end
            
            
        end
    end
    
    
end

%double check moviet0 was properly assigned
assert(~isnan(movieT0));

Frame_Times = Frame_Times - movieT0;

end