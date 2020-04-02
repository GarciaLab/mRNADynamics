
function Particles = addPositionsToParticles(Particles, Spots, currentChannel)


    %AR: add x, y and z positions to the Particles structure. This was
    %originally only added by addParticlePosition, but it it's more useful
    %early on in the pipeline.
    
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