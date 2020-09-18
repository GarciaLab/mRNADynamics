function [CurrentNucleus, CurrentFrame, ManualZFlag] =...
    goEarliestUncheckedNucleus(CurrentNucleus, schnitzcells)
%GONEXTPARTICLE Summary of this function goes here
%   Detailed explanation goes here

%lineFit = 0; % the initial rise was not fitted!
%fitApproved = 0; % the initial rise fit was not approved!
numNuclei = length(schnitzcells);
ApprovedNuclei = [schnitzcells.Approved];
CheckedNuclei = [schnitzcells.Checked];
nuclei_idx = 1:numNuclei;
idx_to_check = nuclei_idx(~CheckedNuclei);
NextNucleus_idx=find(idx_to_check ~= CurrentNucleus, 1);
NextNucleus = idx_to_check(NextNucleus_idx);
if isempty(NextNucleus)
    NextNucleus=numNuclei;
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

