function [SpotFilter,Particles]=...
    TransferParticle(Spots,SpotFilter,Particles,CurrentFrame,IndexOutput)

%First, approve the particle
SpotFilter(CurrentFrame,IndexOutput)=1;

 %HG + AR: We deleted this because we don't assign xPos and yPos for the
 %Particles that are detected in the first place. %AR 8/16/2018- reinstated
 %this because it was causing bugs. 
% %Get the position of the this particle
[x,y,z]=SpotsXYZ(Spots(CurrentFrame));
x=x(IndexOutput);
y=y(IndexOutput);
z=z(IndexOutput);

%Add this spot as a new particle to the end of the Particles structure
Particles(end+1).Frame=CurrentFrame;
Particles(end).Index=IndexOutput;
Particles(end).Approved=false;
 Particles(end).FrameApproved=true;

 %HG + AR: We deleted this because we don't assign xPos and yPos for the
 %Particles that are detected in the first place.%AR 8/16/2018- reinstated
 %this because it was causing bugs. 
Particles(end).xPos=x;
Particles(end).yPos=y;
Particles(end).zPos=z;