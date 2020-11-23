function [schnitzcells, CurrentNucleus] = orderNuclei(numNuclei, schnitzcells, CurrentNucleus, ReorderOrientation)
%ORDERPARTICLES Summary of this function goes here
%   Detailed explanation goes here

%Order particles by the earliest frame they appear at. This makes the
%tracking a lot easier!
clear schnitz_params
schnitz_params = zeros(numNuclei, 3, 'double');
for i=1:numNuclei
    schnitz_params(i, 1)=schnitzcells(i).frames(1);
    schnitz_params(i, 2) = schnitzcells(i).cenx(1);
    schnitz_params(i, 3) = schnitzcells(i).ceny(1);
end
if ReorderOrientation == 0
    schnitz_params(:,2) = -1*schnitz_params(:,2);
end
[~,Permutations]=sortrows(schnitz_params);
schnitzcells=schnitzcells(Permutations);
CurrentNucleus = find(Permutations == CurrentNucleus, 1);
end

