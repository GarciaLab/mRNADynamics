
function Particles = addPositionsToParticles(Particles, Spots, currentChannel)


    %add x, y and z positions to the Particles structure.
    
    for particle=1:length(Particles{currentChannel})
        for frame=1:length(Particles{currentChannel}(particle).Frame)
            [x,y,z]=SpotsXYZ(Spots{currentChannel}(Particles{currentChannel}(particle).Frame(frame)));
            if ~isempty(x)
                Particles{currentChannel}(particle).xPos(frame)=...
                    x(Particles{currentChannel}(particle).Index(frame));
                
                Particles{currentChannel}(particle).yPos(frame)=...
                    y(Particles{currentChannel}(particle).Index(frame));
                
                Particles{currentChannel}(particle).zPos(frame)=...
                    z(Particles{currentChannel}(particle).Index(frame));
            end
        end
    end
    


end