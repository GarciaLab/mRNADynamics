function Frame_Times = getFrameTimesFromBioFormats(LIFMeta)

import ome.xml.model.primitives.Timestamp;
import ome.units.UNITS;
import ome.units.quantity.Time;
import loci.common.DateTools;

nSeries = LIFMeta.getImageCount();

Frame_Times = [];

for s = 0:nSeries-1
    
    nPlanes = LIFMeta.getPlaneCount(s);
    
    seriesStamp = LIFMeta.getImageAcquisitionDate(s);
    seriesStamp_char = seriesStamp.getValue();
    t0 = DateTools.getTime(seriesStamp_char, DateTools.ISO8601_FORMAT) / 1000; %unix epoch ms -> s
    
    if s == 0
        t00 = t0;
    end
    
    ft = [];
    
    for p = 0:nPlanes-1
        
        dt = LIFMeta.getPlaneDeltaT(s, p);
        if ~isempty(dt)
            ft = [ft, t0 + double(dt.value)];
        end
        
        
    end
    
    
    Frame_Times = [Frame_Times, ft];
    
end

Frame_Times = Frame_Times - t00;

end