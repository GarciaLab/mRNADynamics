function NucleusOutput = identifyNucleus(schnitzcells, CurrentFrame, ...
    UseHistoneOverlay, ConnectPosition)
%IDENTIFYPARTICLE Summary of this function goes here
%   Detailed explanation goes here

%numNuclei = length(schnitzcells);
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

end

