function [Particles,schnitzcells]=CheckSchnitzLineage(Particles,CurrentParticle,schnitzcells,CurrentFrame,...
    Overlay)


%This function checks that every frame within a particle points at the one
%schnitz. It's useful, for now, when connecting to a fad2 particle.

%First check if the daughter nuclei exist in this frame. We have to do this
%because of the puncta that persist beyond a nucleus existence.
DaughterFrames=[];
if schnitzcells(Particles(CurrentParticle).Nucleus).E~=0
    DaughterFrames=[DaughterFrames,schnitzcells(schnitzcells(Particles(CurrentParticle).Nucleus).E).frames];
end
if schnitzcells(Particles(CurrentParticle).Nucleus).D~=0
    DaughterFrames=[DaughterFrames,schnitzcells(schnitzcells(Particles(CurrentParticle).Nucleus).D).frames];
end


if (~sum((schnitzcells(Particles(CurrentParticle).Nucleus).frames)==CurrentFrame))&...
        ~sum((DaughterFrames-1)==CurrentFrame)
    display('Error with the nuclear tracking. Select corresponding nucleus')
    
    %Ask the user to click on the corresponding nucleus
    figure(Overlay)
    SelectedNucleus=ginput(1);

    %Find which schnitz this corresponds to
    SchnitzSuspect=[];
    xPosSuspect=[];
    yPosSuspect=[];
    for i=1:length(schnitzcells)
        if sum(schnitzcells(i).frames==CurrentFrame)
            SchnitzSuspect=[SchnitzSuspect,i];
            xPosSuspect=[xPosSuspect,...
                schnitzcells(i).cenx(find((schnitzcells(i).frames==CurrentFrame))];
            yPosSuspect=[yPosSuspect,...
                schnitzcells(i).ceny(find((schnitzcells(i).frames)==CurrentFrame))];
        end
    end
    
    
    %Find the closest one to the point where we clicked
    Distance=sqrt((SelectedNucleus(1)-xPosSuspect).^2+(SelectedNucleus(2)-yPosSuspect).^2);
    [MinValue,ClosestNucleusIndex]=min(Distance);
    
    
    %Join the corresponding schnitz
    [Particles,schnitzcells]=JoinSchnitz(Particles,schnitzcells,Particles(CurrentParticle).Nucleus,...
        SchnitzSuspect(ClosestNucleusIndex))
end    