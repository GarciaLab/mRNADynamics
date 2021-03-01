function PlotParticleTrace(cptState, plotTraceSettings, noSpline, UseTwin)
% noSpline displays the particle trace as well as a montage of the images
if ~exist('UseTwin', 'var')
    UseTwin = false;
end
if ~UseTwin
    CurrentParticle = cptState.CurrentParticle;
    Particles = cptState.getCurrentChannelParticles();
    Spots = cptState.getCurrentChannelSpots();
    
    GetParticleTrace(CurrentParticle, Particles, Spots, plotTraceSettings, noSpline);
    
    if ~plotTraceSettings.UseCompiledParticles
        cptState.Frames = Particles(CurrentParticle).Frame;
    else
        cptState.Frames = Particles(CurrentParticle).FlaggingInfo.TrueFrames;
    end
else
    TwinParticle = cptState.TwinParticle;
    Particles = cptState.getCurrentChannelParticles();
    Spots = cptState.getCurrentChannelSpots();
    
    GetParticleTrace(TwinParticle, Particles, Spots, plotTraceSettings, noSpline);
    
    if ~plotTraceSettings.UseCompiledParticles
        cptState.Frames = Particles(TwinParticle).Frame;
    else
        cptState.Frames = Particles(TwinParticle).FlaggingInfo.TrueFrames;
    end
end
%{
    for i=1:length(cptState.Frames)
        
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
%}
end
