function [lineFit, CurrentParticle, CurrentFrame, ManualZFlag, DisplayRange, TwinParticle] =...
    goNextUncheckedParticle(CurrentParticle, CurrentChannelIndex, HideApprovedFlag, Particles, UseTwinTraces)
%GONEXTPARTICLE Summary of this function goes here
%   Detailed explanation goes here

if ~exist('UseTwinTraces', 'var')
    UseTwinTraces = false;
end

lineFit = 0; % the initial rise was not fitted!
fitApproved = 0; % the initial rise fit was not approved!
numParticles = length(Particles{CurrentChannelIndex});


ParticleApprovals = [Particles{CurrentChannelIndex}(:).Approved];
ParticleIndices = 1:numParticles;
hasNucleus = zeros(1, numParticles);
for i= ParticleIndices
    if ~isempty(Particles{CurrentChannelIndex}(i).Nucleus)
        hasNucleus(i) = 1;
    end
end
NextParticle = find(ParticleApprovals == 0 & ParticleIndices > CurrentParticle & hasNucleus == 1, 1);
if isempty(NextParticle)
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


[CurrentParticle,CurrentFrame, ManualZFlag, TwinParticle] = ...
    changeParticle(NextParticle, Particles, numParticles, CurrentChannelIndex, UseTwinTraces);

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

