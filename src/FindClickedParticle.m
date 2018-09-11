function [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,Spots,Particles)

%Find the particle index of the particle the user clicked on



%Find the closest particle and its index
[x,y]=SpotsXYZ(Spots(CurrentFrame));
Distances=sqrt((x-ConnectPosition(1)).^2+(y-ConnectPosition(2)).^2);
[~,IndexOutput]=min(Distances);

%Now, look for the particle in the Particles structure
if isempty(Distances)
    ParticleOutput=[];
    IndexOutput=[];
else
    for i=1:length(Particles)
        FrameFilter=(Particles(i).Frame==CurrentFrame);
        if Particles(i).Index(FrameFilter)==IndexOutput
            ParticleOutput=i;
        end
    end
end