function [lineFit, CurrentParticle, CurrentFrame, ManualZFlag, DisplayRange] =...
    goNextParticle(CurrentParticle, CurrentChannelIndex, HideApprovedFlag, Particles)
%GONEXTPARTICLE Summary of this function goes here
%   Detailed explanation goes here

lineFit = 0; % the initial rise was not fitted!
fitApproved = 0; % the initial rise fit was not approved!
numParticles = length(Particles{CurrentChannelIndex});
NextParticle=CurrentParticle+1;

if NextParticle>numParticles
    NextParticle=numParticles;
end


%Mode 1 - skip approved or flagged traces
while (HideApprovedFlag)==1&&(NextParticle<numParticles)&&...
        ((Particles{CurrentChannelIndex}(NextParticle).Approved==1)||(Particles{CurrentChannelIndex}(NextParticle).Approved==-1)||...
        (Particles{CurrentChannelIndex}(NextParticle).Approved==2))
    NextParticle=NextParticle+1;
end

%Mode 2 - skip approved traces
while ((HideApprovedFlag)==2)&&(NextParticle<numParticles)&&...
        ((Particles{CurrentChannelIndex}(NextParticle).Approved==1)||(Particles{CurrentChannelIndex}(NextParticle).Approved==2))
    NextParticle=NextParticle+1;
end


[CurrentParticle,CurrentFrame, ManualZFlag] = ...
    changeParticle(NextParticle, Particles, numParticles, CurrentChannelIndex);

DisplayRange=[];

msg = Particles{CurrentChannelIndex}(CurrentParticle).Frame(find(diff(Particles{CurrentChannelIndex}(CurrentParticle).Frame)>1));

if ~isempty(msg)
    %             disp('Missing frames:') %AR 12/3/17- Not sure what this
    %             message is trying to say, so I am silencing it for now.
    %             msg
else
    %do nothing
end
end

