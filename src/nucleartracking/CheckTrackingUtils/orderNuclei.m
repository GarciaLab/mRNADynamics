function schnitzcells = orderNuclei(numNuclei, schnitzcells)
%ORDERPARTICLES Summary of this function goes here
%   Detailed explanation goes here

%Order particles by the earliest frame they appear at. This makes the
%tracking a lot easier!
clear FirstFrame
for i=1:numNuclei
    FirstFrame(i)=schnitzcells(i).frames(1);
end
[~,Permutations]=sort(FirstFrame);
schnitzcells=schnitzcells(Permutations);
end

