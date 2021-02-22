function NucleusOutput = SelectNucleus(schnitzcells, Ellipses, CurrentFrame, ConnectPosition)
%IDENTIFYPARTICLE Summary of this function goes here
%   Detailed explanation goes here

%Find the closest particle and its index
x = Ellipses{CurrentFrame}(:,1);
y = Ellipses{CurrentFrame}(:,2);
Distances=sqrt((x-ConnectPosition(1)).^2+(y-ConnectPosition(2)).^2);
[~,IndexOutput]=min(Distances);
nNuclei = length(schnitzcells);

%Now, look for the particle in the Particles structure
if isempty(Distances)
    NucleusOutput=[];
    IndexOutput=[];
else
    NucleusOutput = Ellipses{CurrentFrame}(IndexOutput,9);
end
