function [fad,fad2,Particles]=...
    Integratefad2Particle(fad,fad2,fad2Position,Particles,CurrentFrame)

%Find the closest particle and its index
[x2,y2]=fad2xyzFit(CurrentFrame,fad2, 'addMargin'); 
Distances=sqrt((x2-fad2Position(1)).^2+(y2-fad2Position(2)).^2);

[Dummy,IndexOutput]=min(Distances);


[fad,fad2,Particles]=...
    TransferParticle(fad,CurrentFrame,[],fad2,CurrentFrame,...
    IndexOutput,Particles);

