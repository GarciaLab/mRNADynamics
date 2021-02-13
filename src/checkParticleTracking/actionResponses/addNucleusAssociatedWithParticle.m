function [Particles, CurrentParticle,CurrentFrame] =...
    addNucleusAssociatedWithParticle(Spots, Particles, PreviousParticle,...
    CurrentFrame, CurrentChannelIndex, UseHistoneOverlay, schnitzcells, ConnectPosition)
%TONEARESTPARTICLE Summary of this function goes here
%   Detailed explanation goes here

try
if ~exist('ConnectPosition', 'var')
    [ConnectPositionx,ConnectPositiony]=ginput(1);
    ConnectPosition = [ConnectPositionx,ConnectPositiony];
end

if ~isempty(ConnectPosition)
    
    display(ConnectPosition);
    
    %Find the closest particle
    [NucleusOutput,~]=FindClickedNucleus(ConnectPosition,CurrentFrame,...
        schnitzcells);
    
    disp(['Clicked nucleus: ',num2str(NucleusOutput)]);
    CurrentParticle = PreviousParticle;
    Particles{CurrentChannelIndex}(CurrentParticle).Nucleus = NucleusOutput;
    
end


catch
    disp('Failed to identify clicked particle. Be careful to click the particle.')
    CurrentParticle = PreviousParticle;
end







end