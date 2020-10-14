function [CurrentNucleus, CurrentFrame, ManualZFlag] =...
    goPreviousNCNucleus(CurrentNucleus, schnitzcells)
%GONEXTPARTICLE Summary of this function goes here
%   Detailed explanation goes here

%lineFit = 0; % the initial rise was not fitted!
%fitApproved = 0; % the initial rise fit was not approved!
numNuclei = length(schnitzcells);
nuclei_idx = 1:numNuclei;
current_nc = schnitzcells(CurrentNucleus).cycle;

nuclear_cycles = [schnitzcells.cycle];
prev_cycle_idx = nuclei_idx(nuclear_cycles == current_nc - 1);
if isempty(prev_cycle_idx) 
    NextNucleus = 1;
else
    NextNucleus = prev_cycle_idx(1);
end




[CurrentNucleus,CurrentFrame, ManualZFlag] = ...
    changeNucleus(NextNucleus, schnitzcells, numNuclei);

% msg = schnitzcells(CurrentNucleus).frames(find(diff(schnitzcells(CurrentNucleus).Frame)>1));
% 
% if ~isempty(msg)
%     %             disp('Missing frames:') %AR 12/3/17- Not sure what this
%     %             message is trying to say, so I am silencing it for now.
%     %             msg
% else
%     %do nothing
% end
end

