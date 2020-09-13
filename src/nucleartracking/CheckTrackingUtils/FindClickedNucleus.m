function [NucleusOutput,IndexOutput]=...
    FindClickedNucleus(ConnectPosition,CurrentFrame,schnitzcells, cntState)

%Find the cellno index of the nucleus/schnitz cell the user clicked on


%Find the closest particle and its index
[x,y,~,cellnos]=NucleiXYZ(schnitzcells, CurrentFrame);
Distances=sqrt((x-ConnectPosition(1)).^2+(y-ConnectPosition(2)).^2);
[~,idx]=min(Distances);
IndexOutput = cellnos(idx);
nNuclei = length(schnitzcells);

%Now, look for the particle in the Particles structure
if isempty(Distances)
    NucleusOutput=[];
    IndexOutput=[];
else
    for p = 1 : nNuclei
        nucleus = schnitzcells(p);
        FrameFilter= nucleus.frames == CurrentFrame;
        if schnitzcells(p).cellno(FrameFilter) == IndexOutput
            NucleusOutput = p;
            break;
        end
    end
end