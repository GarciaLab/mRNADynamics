function Particles = addTracksToParticles(Particles, tracks, CurrentFrame)

% This is a subfunction for tracknuclei_computervision 

    for k = 1:length(tracks)
        if k > length(Particles) 

            Particles(k).frames = CurrentFrame;
            Particles(k).xPos = tracks(k).kalmanFilter.State(1);
            Particles(k).yPos = tracks(k).kalmanFilter.State(4);
            Particles(k).zPos = tracks(k).kalmanFilter.State(7); 
            Particles(k).FramesSinceDetection = tracks(k).consecutiveInvisibleCount;
            Particles(k).Approved = 0;
            Particles(k).Index = k;

        else

            Particles(k).frames(end+1) = CurrentFrame;
            Particles(k).xPos(end+1) = tracks(k).kalmanFilter.State(1);
            Particles(k).yPos(end+1) = tracks(k).kalmanFilter.State(4);
            Particles(k).zPos(end+1) = tracks(k).kalmanFilter.State(7);      
            Particles(k).FramesSinceDetection(end+1) = tracks(k).consecutiveInvisibleCount;

        end
    end


end