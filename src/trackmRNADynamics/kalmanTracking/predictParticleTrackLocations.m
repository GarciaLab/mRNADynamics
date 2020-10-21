function particleTracks = predictParticleTrackLocations(particleTracks)

        for k = 1:length(particleTracks)
                      
            % Predict the current location of the track.           
            predict(particleTracks(k).kalmanFilter);                       
            
        end
    end
