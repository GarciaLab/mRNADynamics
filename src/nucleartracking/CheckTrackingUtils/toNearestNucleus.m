function [CurrentNucleus, CurrentFrame, ManualZFlag] = toNearestNucleus(schnitzcells, ...
    CurrentFrame, UseHistoneOverlay, ConnectPosition)
%TONEARESTNUCLEUS Summary of this function goes here
%   Detailed explanation goes here

numNuclei = length(schnitzcells);
opts = {};
if exist('ConnectPosition', 'var')
    opts = {'ConnectPosition'};
end

NucleusOutput = identifyNucleus(schnitzcells, CurrentFrame, ...
    UseHistoneOverlay, opts{:});

if (floor(NucleusOutput)>0)&(NucleusOutput<=numNuclei)
    [CurrentNucleus,CurrentFrame, ManualZFlag] = ...
    changeNucleus(NucleusOutput, schnitzcells, numNuclei);
end

end

