function [CurrentChannel, PreviousChannel, coatChannel, CurrentParticle] =...
    switchChannels(CurrentChannel, CurrentParticle, Particles, ...
    UseHistoneOverlay, NChannels)
%SWITCHCHANNELS Summary of this function goes here
%   Detailed explanation goes here

%NB Armando 3/28/2020-
%I think this is completely broken. Needs a rewrite using
%CurrentChannelIndex and CurrentChannel

%Update the channel number
PreviousChannel=CurrentChannel;
CurrentChannel=CurrentChannel+1;
if CurrentChannel>NChannels
    CurrentChannel=1;
end

%Update the coatChannel
coatChannels = [1, 2]; % JP: this looks weird, but it was this way in the main loop

coatChannel=coatChannels(CurrentChannel);


%Do we have a histone channel? If so, we can find the particle in
%the next channel corresponding to this nucleus.
numParticlesCurrCh = length(Particles{CurrentChannel});
if UseHistoneOverlay
    %If a particle is associated with this same nucleus in the new
    %channel then change to it
    AssignedNucleusPreviousChannel=Particles{PreviousChannel}(CurrentParticle).Nucleus;
    %Now, find its associated particle
    AssignedNucleusNewChannel=[];
    for i=1:numParticlesCurrCh
        if ~isempty(Particles{CurrentChannel}(i).Nucleus)
            AssignedNucleusNewChannel(i)=Particles{CurrentChannel}(i).Nucleus;
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

