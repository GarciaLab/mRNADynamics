function [CurrentParticle,CurrentFrame, ManualZFlag] = ...
    changeParticle(ParticleNum, Particles, numParticles, CurrentChannel)
%CHANGEPARTICLE Summary of this function goes here
%   Detailed explanation goes here
%     lineFit = 0; % the initial rise was not fitted!
%     fitApproved = 0; % the initial rise fit was not approved!
    CurrentParticle = min(max(ParticleNum, 1), numParticles);
    ManualZFlag = 0;
    CurrentFrame = Particles{CurrentChannel}(CurrentParticle).Frame(1);
end

