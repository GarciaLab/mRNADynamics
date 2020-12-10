function [Particles, CurrentParticle,CurrentFrame, ManualZFlag] =...
    changeParticleAssociatedWithNucleus(Spots, Particles, PreviousParticle,...
    CurrentFrame, CurrentChannelIndex, UseHistoneOverlay, schnitzcells, ConnectPosition)
%TONEARESTPARTICLE Summary of this function goes here
%   Detailed explanation goes here


numParticles = length(Particles{CurrentChannelIndex});
opts = {};
if exist('ConnectPosition', 'var')
    opts = {'ConnectPosition'};
end
try
    ParticleOutput = identifyParticle(Spots, Particles, CurrentFrame, ...
        CurrentChannelIndex, UseHistoneOverlay, schnitzcells, opts{:});
    
    if (floor(ParticleOutput)>0)&(ParticleOutput<=numParticles)
        
        [CurrentParticle,CurrentFrame, ManualZFlag] = ...
            changeParticle(ParticleOutput, Particles, numParticles, CurrentChannelIndex);
    end
catch
    disp('Failed to identify clicked particle. Be careful to click the particle.')
end



Particles{CurrentChannelIndex}(CurrentParticle).Nucleus = Particles{CurrentChannelIndex}(PreviousParticle).Nucleus;
Particles{CurrentChannelIndex}(PreviousParticle).Nucleus = [];
Particles{CurrentChannelIndex}(PreviousParticle).Approved = 0;



end