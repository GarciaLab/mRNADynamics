function [CurrentParticle,CurrentFrame, ManualZFlag, TwinParticle] = ...
    changeParticle(ParticleNum, Particles, numParticles, CurrentChannel, UseTwinTraces)

if ~exist('UseTwinTraces', 'var')
    UseTwinTraces = false;
end


CurrentParticle = min(max(ParticleNum, 1), numParticles);
ManualZFlag = 0;
CurrentFrame = Particles{CurrentChannel}(CurrentParticle).Frame(1);

if UseTwinTraces
    TwinParticle = findTwinParticle(CurrentParticle, Particles, CurrentChannel);
else
    TwinParticle = [];
end
 
end

