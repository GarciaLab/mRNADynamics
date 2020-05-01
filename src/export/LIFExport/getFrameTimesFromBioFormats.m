function Frame_Times = getFrameTimesFromBioFormats(LIFMeta)

import loci.common.DateTools;

nSeries = LIFMeta.getImageCount();

Frame_Times = [];

for series = 0:nSeries-1
    
    nPlanes = LIFMeta.getPlaneCount(series);
    
    seriesStamp = LIFMeta.getImageAcquisitionDate(series);
    seriesStamp_char = seriesStamp.getValue();
    seriesT0 = DateTools.getTime(seriesStamp_char, DateTools.ISO8601_FORMAT) / 1000; %unix epoch ms -> s
    
    if series == 0
        movieT0 = seriesT0;
    end
    
    
    for plane = 0:nPlanes-1
        
        deltaT = LIFMeta.getPlaneDeltaT(series, plane);
        if ~isempty(deltaT)
            Frame_Times = [Frame_Times, seriesT0 + double(deltaT.value)];
        end
        
        
    end
    
        
end

Frame_Times = Frame_Times - movieT0;

end