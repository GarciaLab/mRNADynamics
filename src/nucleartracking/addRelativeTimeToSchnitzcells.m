function [schnitzcells, ncTimes] = addRelativeTimeToSchnitzcells(schnitzcells, FrameInfo, ncFrames)


    ncFrames(ncFrames==0) = 1;
    ind = find(isnan(ncFrames));
    for k = 1:length(ind)
        ncFrames(ind(k)) = ncFrames(ind(k)-1); %if the last nc was 'NaN' change it to the previous nc frame
    end
    time = [FrameInfo.Time]/60; %frame times in minutes 
    ncTimes = time(ncFrames);
    
    for s = 1:length(schnitzcells)
        schnitzcells(s).timeSinceAnaphase = time(schnitzcells(s).frames) - ncTimes(schnitzcells(s).cycle);
        % if the movie didn't capture the anaphase of this schnitz' cycle,
        % we want to keep track of that
        if ncTimes(schnitzcells(s).cycle) == 0 
            schnitzcells(s).catchedAnaphase = 0;
        else
            schnitzcells(s).catchedAnaphase = 1;
        end
    end
    
end