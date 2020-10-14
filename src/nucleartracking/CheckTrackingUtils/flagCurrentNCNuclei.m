function [CurrentNucleus, CurrentFrame, schnitzcells] =...
    flagCurrentNCNuclei(CurrentNucleus, schnitzcells)
%GONEXTPARTICLE Summary of this function goes here
%   Detailed explanation goes here

%lineFit = 0; % the initial rise was not fitted!
%fitApproved = 0; % the initial rise fit was not approved!
numNuclei = length(schnitzcells);

nuclei_idx = 1:numNuclei;
current_nc = schnitzcells(CurrentNucleus).cycle;
CheckedStatus = [schnitzcells(:).Checked];
nuclear_cycles = [schnitzcells.cycle];
this_cycle_idx = nuclei_idx((nuclear_cycles == current_nc) & ~Checked);

for i=1:length(this_cycle_idx)
    idx = this_cycle_idx(i);
    schnitzcells(idx).Checked = 1;
    schnitzcells(idx).Approved = 0;
    schnitzcells(idx).Flag = 7;
end
next_cycle_idx = nuclei_idx(nuclear_cycles == current_nc + 1);
if isempty(next_cycle_idx) 
    NextNucleus = numNuclei;
else
    NextNucleus = next_cycle_idx(1);
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

