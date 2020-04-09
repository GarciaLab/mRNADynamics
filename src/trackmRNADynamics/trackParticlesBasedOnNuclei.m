function [Particles, SpotFilter] = trackParticlesBasedOnNuclei(...
    ...
    PreProcPath, Prefix, CurrentFrame, NDigits, app, nucAxes, Ellipses, ...
    ExperimentType, Channel, schnitzcells,...
    Particles, Spots, SpotFilter, PixelSize,...
    SearchRadius, retrack, displayFigures, thisExperiment)


% 
% trackParticlesBasedOnNuclei(PreProcPath, Prefix,...
%                 CurrentFrame, NDigits, app, nucAxes, Ellipses, ...
%                 ExperimentType, Channel, scurrentChannelitzcells, Particles, Spots,...
%                 SpotFilter, PixelSize, SearchRadius, retrack, displayFigures, thisExperiment)
            
            
            


if displayFigures
    hisImage = openHistoneImage(Prefix, PreProcPath, CurrentFrame, NDigits);
    
    if ~isempty(app)
        ax2 = app{2};
    else
        ax2 = nucAxes;
    end
    
    imshow(hisImage, [], 'Border', 'Tight', 'Parent', ax2, 'InitialMagnification', 'fit')
    hold(ax2, 'on')
    PlotHandle = [];
    [NEllipses, ~] = size(Ellipses{CurrentFrame});
    
    for EllipsesIndex = 1:NEllipses
        PlotHandle = [PlotHandle, ellipse(...
            Ellipses{CurrentFrame}(EllipsesIndex, 3), ...
            Ellipses{CurrentFrame}(EllipsesIndex, 4), ...
            Ellipses{CurrentFrame}(EllipsesIndex, 5), ...
            Ellipses{CurrentFrame}(EllipsesIndex, 1) + 1, ...
            Ellipses{CurrentFrame}(EllipsesIndex, 2) + 1, ...
            [], [], ax2)];
        
        text(ax2, Ellipses{CurrentFrame}(EllipsesIndex, 1) + 1,...
            Ellipses{CurrentFrame}(EllipsesIndex, 2) + 1, ...
            num2str(EllipsesIndex), 'BackgroundColor', [.7 .9 .7]);
    end
    
    set(PlotHandle, 'Color', 'r')
    hold(ax2, 'off')
    title(ax2, CurrentFrame)
    drawnow
end

if strcmp(ExperimentType, '2spot')
    SpotsPerNucleus = 2;
else
    SpotsPerNucleus = 1;
end

[Particles{Channel}, SpotFilter{Channel}] = AssignParticle2Nucleus(...
    schnitzcells, Ellipses, ...
    Particles{Channel}, Spots{Channel}, SpotFilter{Channel},...
    CurrentFrame, PixelSize, SpotsPerNucleus, retrack);

end