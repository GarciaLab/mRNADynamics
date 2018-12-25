function [lineFit, CurrentParticle, CurrentFrame, ManualZFlag, DisplayRange] =...
    goNextParticle(CurrentParticle, CurrentChannel, HideApprovedFlag, Particles)
%GONEXTPARTICLE Summary of this function goes here
%   Detailed explanation goes here

lineFit = 0; % the initial rise was not fitted!
numParticles = length(Particles{CurrentChannel});
NextParticle=CurrentParticle+1;

if NextParticle>numParticles
    NextParticle=numParticles;
end


%Mode 1 - skip approved or flagged traces
while (HideApprovedFlag)==1&&(NextParticle<numParticles)&&...
        ((Particles{CurrentChannel}(NextParticle).Approved==1)||(Particles{CurrentChannel}(NextParticle).Approved==-1)||...
        (Particles{CurrentChannel}(NextParticle).Approved==2))
    NextParticle=NextParticle+1;
end

%Mode 2 - skip approved traces
while ((HideApprovedFlag)==2)&&(NextParticle<numParticles)&&...
        ((Particles{CurrentChannel}(NextParticle).Approved==1)||(Particles{CurrentChannel}(NextParticle).Approved==2))
    NextParticle=NextParticle+1;
end


[CurrentParticle,CurrentFrame, ManualZFlag] = ...
    changeParticle(NextParticle, Particles, numParticles, CurrentChannel);

DisplayRange=[];

msg = Particles{CurrentChannel}(CurrentParticle).Frame(find(diff(Particles{CurrentChannel}(CurrentParticle).Frame)>1));

if ~isempty(msg)
    %             disp('Missing frames:') %AR 12/3/17- Not sure what this
    %             message is trying to say, so I am silencing it for now.
    %             msg
else
    %do nothing
end
end

