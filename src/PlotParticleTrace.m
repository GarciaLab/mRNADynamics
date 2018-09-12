function [Frames,AmpIntegral,AmpGaussian,AmpIntegral3,AmpIntegral5,...
    ErrorIntegral, ErrorIntegral3, ErrorIntegral5,backGround3]=...
    PlotParticleTrace(CurrentParticle,Particles,Spots)

%This displays the particle trace as well as a montage of the images

%V2: Removed the ImageSnippet for now. I was running into the borders of
%the image.

[Frame,AmpIntegral,AmpIntegral3,AmpIntegral5,AmpGaussian,Offset,...
    ErrorIntegral,ErrorGauss,optFit,FitType,ErrorIntegral3, ErrorIntegral5,backGround3]=GetParticleTrace(CurrentParticle,Particles,Spots);
Frames=Particles(CurrentParticle).Frame;
Indexes=Particles(CurrentParticle).Index;



for i=1:length(Frames)
    CurrentFrame=Frames(i);

    %Get the coordinates taking the margins into account
    [x,y,z]=SpotsXYZ(Spots(CurrentFrame));
    
    %Pull out the right particle if it exists in this frame
    CurrentParticleIndex=Particles(CurrentParticle).Index(Particles(CurrentParticle).Frame==CurrentFrame);
    
    %This is the position of the current particle
    xTrace=round(x(CurrentParticleIndex));
    yTrace=round(y(CurrentParticleIndex));
    zTrace=round(z(CurrentParticleIndex));
end

