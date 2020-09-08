function [NucleusOutput,IndexOutput]=...
    FindClickedNucleus(ConnectPosition,CurrentFrame,schnitzcells)

%Find the particle index of the particle the user clicked on


%Find the closest particle and its index
[x,y]=NucleiXYZ(schnitzcells, CurrentFrame);
Distances=sqrt((x-ConnectPosition(1)).^2+(y-ConnectPosition(2)).^2);
[~,IndexOutput]=min(Distances);

nNuclei = length(schnitzcxells);

%Now, look for the particle in the Particles structure
if isempty(Distances)
    NucleusOutput=[];
    IndexOutput=[];
else
    for p = 1 : nNuclei
        nucleus = schnitzcells(p);
        FrameFilter= nucleus.frames == CurrentFrame;
        if schnitzcells.cellno(FrameFilter) == IndexOutput
            NucleusOutput = p;
            break;
        end
    end
end