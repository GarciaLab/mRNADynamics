function [CurrentParticle, CurrentFrame, ManualZFlag] = toNearestParticle(Spots, ...
    Particles, CurrentFrame, CurrentChannel, UseHistoneOverlay, schnitzcells, ConnectPosition)
%TONEARESTPARTICLE Summary of this function goes here
%   Detailed explanation goes here

numParticles = length(Particles{CurrentChannel});
opts = {};
if exist('ConnectPosition', 'var')
    opts = {'ConnectPosition'};
end

ParticleOutput = identifyParticle(Spots, Particles, CurrentFrame, ...
    CurrentChannel, UseHistoneOverlay, schnitzcells, opts{:});

if (floor(ParticleOutput)>0)&(ParticleOutput<=numParticles)
    [CurrentParticle,CurrentFrame, ManualZFlag] = ...
    changeParticle(ParticleOutput, Particles, numParticles, CurrentChannel);
end

end

