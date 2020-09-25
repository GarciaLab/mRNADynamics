function [CurrentNucleus,CurrentFrame, ManualZFlag] = ...
    changeNucleus(NucleusNum, schnitzcells, numNuclei)

    CurrentNucleus = min(max(NucleusNum, 1), numNuclei);
    ManualZFlag = 0;
    
    CurrentFrame = schnitzcells(CurrentNucleus).frames(1);


end

