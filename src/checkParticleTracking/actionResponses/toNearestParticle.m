function [CurrentParticle, CurrentFrame, ManualZFlag] = toNearestParticle(Spots, ...
    Particles, CurrentFrame, CurrentChannel, UseHistoneOverlay, schnitzcells)
%TONEARESTPARTICLE Summary of this function goes here
%   Detailed explanation goes here

numParticles = length(Particles{CurrentChannel});

ParticleOutput = identifyParticle(Spots, Particles, CurrentFrame, ...
    CurrentChannel, UseHistoneOverlay, schnitzcells);

if (floor(ParticleOutput)>0)&&(ParticleOutput<=numParticles)
    [CurrentParticle,CurrentFrame, ManualZFlag] = ...
    changeParticle(ParticleOutput, Particles, numParticles, CurrentChannel);
end

end

