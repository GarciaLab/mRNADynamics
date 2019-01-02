function Particles = orderParticles(numParticles, CurrentChannel, Particles)
%ORDERPARTICLES Summary of this function goes here
%   Detailed explanation goes here

%Order particles by the earliest frame they appear at. This makes the
%tracking a lot easier!
clear FirstFrame
for i=1:numParticles
    FirstFrame(i)=Particles{CurrentChannel}(i).Frame(1);
end
[~,Permutations]=sort(FirstFrame);
Particles{CurrentChannel}=Particles{CurrentChannel}(Permutations);
end

