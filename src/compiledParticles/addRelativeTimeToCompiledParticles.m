function CompiledParticles = addRelativeTimeToCompiledParticles(...
    schnitzcells, CompiledParticles, ncFrames, FrameInfo)

    ncFrames(ncFrames==0) = 1;
    ind = find(isnan(ncFrames));
    ncFrames(ind) = ncFrames(ind-1);
    time = [FrameInfo.Time]/60; %frame times in minutes 
    ncTimes = time(ncFrames);

for p = 1:length(CompiledParticles)
    
    s = CompiledParticles(p).Nucleus;
    
    CompiledParticles(p).timeSinceAnaphase =...
        time(CompiledParticles(p).Frame) - ncTimes(schnitzcells(s).cycle);
    
    
    %why did this particle start at the beginning of anaphase? it's 
    %suspicious.
     if CompiledParticles(p).timeSinceAnaphase == 0 & schnitzcells(s).cycle == 12
        warning('stop')
     end
    
end


end