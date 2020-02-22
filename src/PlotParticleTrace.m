function PlotParticleTrace(cptState, noSpline)
    % noSpline displays the particle trace as well as a montage of the images

    CurrentParticle = cptState.CurrentParticle;
    Particles = cptState.getCurrentChannelParticles();
    Spots = cptState.getCurrentChannelSpots();

    [Frame,AmpIntegral,AmpIntegral3,AmpGaussian,Offset,...
        ErrorIntegral,ErrorGauss,~, ~,ErrorIntegral3,...
        backGround3, AmpIntegralGauss3D, ErrorIntegralGauss3D]=...
        ...
        GetParticleTrace(CurrentParticle, Particles, Spots, noSpline);

    cptState.Frames = Particles(CurrentParticle).Frame;
    Indexes = Particles(CurrentParticle).Index;

    for i=1:length(Frames)
        
        CurrentFrame = Frames(i);

        %Get the coordinates taking the margins into account
        [x,y,z] = SpotsXYZ(Spots(CurrentFrame));
        
        %Pull out the right particle if it exists in this frame
        CurrentParticleIndex = Particles(CurrentParticle).Index(Particles(CurrentParticle).Frame == CurrentFrame);
        
        %This is the position of the current particle
        xTrace = round(x(CurrentParticleIndex));
        yTrace = round(y(CurrentParticleIndex));
        zTrace = round(z(CurrentParticleIndex));
    end
end
