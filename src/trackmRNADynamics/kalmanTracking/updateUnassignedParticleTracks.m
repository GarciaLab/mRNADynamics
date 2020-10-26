    function particleTracks = updateUnassignedParticleTracks(particleTracks, ...
                    unassignedTracks, maxUnobservedFrames)  
                  
        for i = 1:length(unassignedTracks)
            ind = unassignedTracks(i);
            particleTracks(ind).age = particleTracks(ind).age + 1;
            particleTracks(ind).consecutiveInvisibleCount = ...
                particleTracks(ind).consecutiveInvisibleCount + 1;
            particleTracks(ind).cappedFlag = particleTracks(ind).consecutiveInvisibleCount > maxUnobservedFrames;
        end
    end