function PlotParticleTrace(cptState, plotTraceSettings, noSpline)
    % noSpline displays the particle trace as well as a montage of the images

    CurrentParticle = cptState.CurrentParticle;
    Particles = cptState.getCurrentChannelParticles();
    Spots = cptState.getCurrentChannelSpots();

    GetParticleTrace(CurrentParticle, Particles, Spots, plotTraceSettings, noSpline);
    
    if ~plotTraceSettings.UseCompiledParticles
        cptState.Frames = Particles(CurrentParticle).Frame;
    else
        cptState.Frames = Particles(CurrentParticle).FlaggingInfo.TrueFrames;
    end
%{    
    for i=1:length(cptState.Frames)
        
        CurrentFrame = Frames(i);

        %Get the coordinates taking the margins into account
        [x,y,z] = getSpotsXYZ(Spots(CurrentFrame));
        
        %Pull out the right particle if it exists in this frame
        CurrentParticleIndex = Particles(CurrentParticle).Index(Particles(CurrentParticle).Frame == CurrentFrame);
        
        %This is the position of the current particle
        xTrace = round(x(CurrentParticleIndex));
        yTrace = round(y(CurrentParticleIndex));
        zTrace = round(z(CurrentParticleIndex));
    end
%}
end
