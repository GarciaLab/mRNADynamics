    function particleTracks = updateUnassignedParticleTracks(particleTracks, ...
                    unassignedTracks, maxUnobservedFrames, trackingOptions, CurrentFrame,fwd_flag)  
                  
        for i = 1:length(unassignedTracks)
            ind = unassignedTracks(i);
            particleTracks(ind).age = particleTracks(ind).age + 1;
            particleTracks(ind).consecutiveInvisibleCount = ...
                particleTracks(ind).consecutiveInvisibleCount + 1;
            % decide whether or not to cap the particle
            
            if fwd_flag
              cap1_flag = CurrentFrame-particleTracks(ind).lastFrame-1 > maxUnobservedFrames;
              cap2_flag = trackingOptions.ncVec(particleTracks(ind).lastFrame) ~= trackingOptions.ncVec(CurrentFrame);
            else
              cap1_flag = particleTracks(ind).firstFrame-CurrentFrame-1 > maxUnobservedFrames;
              cap2_flag = trackingOptions.ncVec(particleTracks(ind).firstFrame) ~= trackingOptions.ncVec(CurrentFrame);
            end
%             if cap1_flag && CurrentFrame > 150
%                 error('wtf')
%             end
            particleTracks(ind).cappedFlag = cap1_flag || cap2_flag;
        end
    end