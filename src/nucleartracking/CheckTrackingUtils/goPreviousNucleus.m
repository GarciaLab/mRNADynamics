function [CurrentNucleus, CurrentFrame, ManualZFlag, DisplayRange] =...
    goPreviousNucleus(CurrentNucleus,  HideApprovedFlag, schnitzcells)
%GOPREVIOUSPARTICLE Summary of this function goes here
%   Detailed explanation goes here

%lineFit = 0; % the initial rise was not fitted!
%fitApproved = 0; % the initial rise fit was not approved!
numNuclei = length(schnitzcells);
NextNucleus=CurrentNucleus-1;

%Mode 1 - show non-flagged traces
while (HideApprovedFlag)==1 && (NextNucleus>1) &&...
        ((schnitzcells(NextNucleus).Approved==1) || (schnitzcells(NextNucleus).Approved==-1) ||...
        (schnitzcells(NextNucleus).Approved==2))
    NextNucleus=NextNucleus-1;
    if NextNucleus<1
        NextNucleus=1;
    end
end


%Mode 2 - show disapproved traces
while ((HideApprovedFlag)==2)&&(NextNucleus>1)&&...
        ((schnitzcells(NextNucleus).Approved==1)||(schnitzcells(NextNucleus).Approved==2))
    NextNucleus=NextNucleus-1;
    if NextNucleus<1
        NextNucleus=1;
    end
end


if NextNucleus<1
    NextNucleus=CurrentNucleus;
end


[CurrentNucleus,CurrentFrame, ManualZFlag] = ...
    changeNucleus(NextNucleus, schnitzcells, numNuclei);


DisplayRange=[];
end

