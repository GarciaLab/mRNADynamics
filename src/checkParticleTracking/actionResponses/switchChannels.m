function [CurrentChannel, PreviousChannel, CurrentParticle,...
    CurrentChannelIndex, PreviousChannelIndex] =...
    switchChannels(CurrentChannel, CurrentParticle, Particles, ...
    UseHistoneOverlay, NChannels, CurrentChannelIndex, PreviousChannelIndex)
%SWITCHCHANNELS Summary of this function goes here
%   Detailed explanation goes here

%NB Armando 3/28/2020-
%I think this is completely broken. Needs a rewrite using
%CurrentChannelIndex and CurrentChannel

coatChannel = nan;

%Update the channel number
PreviousChannel=CurrentChannel;
PreviousChannelIndex = CurrentChannelIndex;
CurrentChannel=CurrentChannel+1;
CurrentChannelIndex = CurrentChannelIndex + 1;
if CurrentChannel>NChannels
    CurrentChannel=1;
    CurrentChannelIndex = 1;
end


%Do we have a histone channel? If so, we can find the particle in
%the next channel corresponding to this nucleus.
numParticlesCurrCh = length(Particles{CurrentChannelIndex});
if UseHistoneOverlay
    %If a particle is associated with this same nucleus in the new
    %channel then change to it
    AssignedNucleusPreviousChannel=Particles{PreviousChannelIndex}(CurrentParticle).Nucleus;
    %Now, find its associated particle
    AssignedNucleusNewChannel=[];
    for i=1:numParticlesCurrCh
        if ~isempty(Particles{CurrentChannelIndex}(i).Nucleus)
            AssignedNucleusNewChannel(i)=Particles{CurrentChannelIndex}(i).Nucleus;
        else
            AssignedNucleusNewChannel(i)=nan;
        end
    end

    if ~isempty(find(AssignedNucleusNewChannel==AssignedNucleusPreviousChannel))
        CurrentParticle=find(AssignedNucleusNewChannel==AssignedNucleusPreviousChannel);
    end

    %If we don't have a histone channel, go for the same particle
    %number in the new channel or for the last particle
elseif numParticlesCurrCh<CurrentParticle
    CurrentParticle=numParticlesCurrCh;
end

end

