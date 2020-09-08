function [CurrentNucleus, CurrentFrame, ManualZFlag, DisplayRange] =...
    goNextNucleus(CurrentNucleus,HideApprovedFlag, schnitzcells)
%GONEXTPARTICLE Summary of this function goes here
%   Detailed explanation goes here

%lineFit = 0; % the initial rise was not fitted!
%fitApproved = 0; % the initial rise fit was not approved!
numNuclei = length(schnitzcells);
NextNucleus=CurrentNucleus+1;

if NextNucleus>numNuclei
    NextNucleus=numNuclei;
end


%Mode 1 - skip approved or flagged traces
while (HideApprovedFlag)==1&&(NextNucleus<numNuclei)&&...
        ((schnitzcells(NextNucleus).Approved==1)||(schnitzcells(NextNucleus).Approved==-1)||...
        (schnitzcells(NextNucleus).Approved==2))
    NextNucleus=NextNucleus+1;
end

%Mode 2 - skip approved traces
while ((HideApprovedFlag)==2)&&(NextNucleus<numNuclei)&&...
        ((schnitzcells(NextNucleus).Approved==1)||(schnitzcells(NextNucleus).Approved==2))
    NextNucleus=NextNucleus+1;
end


[CurrentNucleus,CurrentFrame, ManualZFlag] = ...
    changeNucleus(NucleusNum, schnitzcells, numNuclei);
DisplayRange=[];

msg = schnitzcells(CurrentNucleus).frames(find(diff(schnitzcells(CurrentParticle).Frame)>1));

if ~isempty(msg)
    %             disp('Missing frames:') %AR 12/3/17- Not sure what this
    %             message is trying to say, so I am silencing it for now.
    %             msg
else
    %do nothing
end
end

