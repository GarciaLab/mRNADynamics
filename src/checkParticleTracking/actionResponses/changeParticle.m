function [CurrentParticle,CurrentFrame, ManualZFlag] = ...
    changeParticle(ParticleNum, Particles, numParticles, CurrentChannel)

    CurrentParticle = min(max(ParticleNum, 1), numParticles);
    ManualZFlag = 0;
    CurrentFrame = Particles{CurrentChannel}(CurrentParticle).Frame(1);
end

