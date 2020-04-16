function [schnitzcells, ncTimes] = addRelativeTimeToSchnitzcells(schnitzcells, FrameInfo, ncFrames)


    ncFrames(ncFrames==0) = 1;
    ind = find(isnan(ncFrames));
    ncFrames(ind) = ncFrames(ind-1);
    time = [FrameInfo.Time]/60; %frame times in minutes 
    ncTimes = time(ncFrames);
    
    for s = 1:length(schnitzcells)
        schnitzcells(s).timeSinceAnaphase = time(schnitzcells(s).frames) - ncTimes(schnitzcells(s).cycle);
    end
    
end