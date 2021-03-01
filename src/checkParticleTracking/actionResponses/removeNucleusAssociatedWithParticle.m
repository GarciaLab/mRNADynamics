function [Particles, CurrentParticle,CurrentFrame] =...
    removeNucleusAssociatedWithParticle(Spots, Particles, CurrentParticle,...
    CurrentFrame, CurrentChannelIndex, UseHistoneOverlay, schnitzcells, ConnectPosition)
%TONEARESTPARTICLE Summary of this function goes here
%   Detailed explanation goes here


Particles{CurrentChannelIndex}(CurrentParticle).Nucleus = [];
Particles{CurrentChannelIndex}(CurrentParticle).Schnitz = NaN;







end