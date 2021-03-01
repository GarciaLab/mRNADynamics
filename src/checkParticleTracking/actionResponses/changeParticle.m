function [CurrentParticle,CurrentFrame, ManualZFlag, TwinParticle] = ...
    changeParticle(ParticleNum, Particles, numParticles, CurrentChannel)

    CurrentParticle = min(max(ParticleNum, 1), numParticles);
    ManualZFlag = 0;
    CurrentFrame = Particles{CurrentChannel}(CurrentParticle).Frame(1);
    
    TwinParticle = findTwinParticle(CurrentParticle, Particles, CurrentChannel);
end

