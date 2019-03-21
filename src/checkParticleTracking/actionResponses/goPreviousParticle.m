function [lineFit, CurrentParticle, CurrentFrame, ManualZFlag, DisplayRange] =...
    goPreviousParticle(CurrentParticle, CurrentChannel, HideApprovedFlag, Particles)
%GOPREVIOUSPARTICLE Summary of this function goes here
%   Detailed explanation goes here

lineFit = 0; % the initial rise was not fitted!
fitApproved = 0; % the initial rise fit was not approved!
numParticles = length(Particles{CurrentChannel});
NextParticle=CurrentParticle-1;

%Mode 1 - show non-flagged traces
while (HideApprovedFlag)==1 && (NextParticle>1) &&...
        ((Particles{CurrentChannel}(NextParticle).Approved==1) || (Particles{CurrentChannel}(NextParticle).Approved==-1) ||...
        (Particles{CurrentChannel}(NextParticle).Approved==2))
    NextParticle=NextParticle-1;
    if NextParticle<1
        NextParticle=1;
    end
end


%Mode 2 - show disapproved traces
while ((HideApprovedFlag)==2)&&(NextParticle>1)&&...
        ((Particles{CurrentChannel}(NextParticle).Approved==1)||(Particles{CurrentChannel}(NextParticle).Approved==2))
    NextParticle=NextParticle-1;
    if NextParticle<1
        NextParticle=1;
    end
end


if NextParticle<1
    NextParticle=CurrentParticle;
end

[CurrentParticle,CurrentFrame, ManualZFlag] = ...
    changeParticle(NextParticle, Particles, numParticles, CurrentChannel);


DisplayRange=[];
end

