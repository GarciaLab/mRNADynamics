function Particles =...
    ...
    addNucleusToParticle(...
    ...
    Particles, CurrentFrame, ...
    CurrentChannel, UseHistoneOverlay, schnitzcells, CurrentParticle)

%addNucleusToParticle Summary of this function goes here
%   Detailed explanation goes here


if ~isempty(schnitzcells) & isfield(Particles{CurrentChannel}, 'xPos')
    
    frames = Particles{CurrentChannel}(CurrentParticle).Frame;
    xPos = Particles{CurrentChannel}(CurrentParticle).xPos(frames == CurrentFrame);
    yPos = Particles{CurrentChannel}(CurrentParticle).yPos(frames == CurrentFrame);
    
    NewNuclei=[xPos, yPos];
    
    %Find which schnitz this corresponds to
    SchnitzSuspect=[];
    xPosSuspect=[];
    yPosSuspect=[];
    for j=1:length(schnitzcells)
        if sum(schnitzcells(j).frames==CurrentFrame)
            SchnitzSuspect=[SchnitzSuspect,j];
            xPosSuspect=[xPosSuspect,...
                schnitzcells(j).cenx((schnitzcells(j).frames)==CurrentFrame)];
            yPosSuspect=[yPosSuspect,...
                schnitzcells(j).ceny((schnitzcells(j).frames)==CurrentFrame)];
        end
    end
    
    %Find the closest one to the point where we clicked
    Distance=sqrt((NewNuclei(1)-xPosSuspect).^2+(NewNuclei(2)-yPosSuspect).^2);
    [~,ClosestNucleusIndex]=min(Distance);
    
    ClickedSchnitz=SchnitzSuspect(ClosestNucleusIndex);
    
    Particles{CurrentChannel}(CurrentParticle).Nucleus = ClickedSchnitz;
else
    disp('Failed to connect particle to a nucleus. Either histone channel isn''t present or AddParticlePosition has not yet been run.');
end

end


