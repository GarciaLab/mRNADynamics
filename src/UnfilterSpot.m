function [SpotFilter,Particles]=...
    UnfilterSpot(Spots,SpotFilter,...
    ClickedSpot,Particles,CurrentFrame)

%Find the closest particle and its index
[x2,y2,z2]=getSpotsXYZ(Spots(CurrentFrame));
Distances=sqrt((x2-ClickedSpot(1)).^2+(y2-ClickedSpot(2)).^2);
[Dummy,IndexOutput]=min(Distances);


[SpotFilter,Particles]=...
    TransferParticle(Spots,SpotFilter,Particles,CurrentFrame,IndexOutput);

