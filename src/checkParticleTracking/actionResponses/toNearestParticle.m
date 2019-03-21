function [CurrentParticle, CurrentFrame, ManualZFlag] = toNearestParticle(Spots, ...
    Particles, CurrentFrame, CurrentChannel, UseHistoneOverlay, schnitzcells, ConnectPosition)
%TONEARESTPARTICLE Summary of this function goes here
%   Detailed explanation goes here

numParticles = length(Particles{CurrentChannel});
if exist('ConnectPosition', 'var')
    ParticleOutput = identifyParticle(Spots, Particles, CurrentFrame, ...
        CurrentChannel, UseHistoneOverlay, schnitzcells, ConnectPosition);
else
    ParticleOutput = identifyParticle(Spots, Particles, CurrentFrame, ...
        CurrentChannel, UseHistoneOverlay, schnitzcells);
end

if (floor(ParticleOutput)>0)&(ParticleOutput<=numParticles)
    [CurrentParticle,CurrentFrame, ManualZFlag] = ...
    changeParticle(ParticleOutput, Particles, numParticles, CurrentChannel);
end

end

