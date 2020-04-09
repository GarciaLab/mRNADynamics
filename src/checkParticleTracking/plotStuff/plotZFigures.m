function MaxZProfile = plotZFigures(zProfileFigAxes, zTraceAxes, ExperimentType, xTrace, cptState, plotTraceSettings, fish, MaxZProfile)

if ~isempty(xTrace)
    % Get the z-DoG profile
    % check to see if Spots contains flag indicating type of integration used
    
    g = [-1 0 1];
    gaussFilter = exp(-g .^ 2 / (2 ));
    zprofinit = zeros(1, cptState.ZSlices);
    
    currentParticleFit = cptState.getCurrentParticleFit();
    
    zprofinit(currentParticleFit.z) = currentParticleFit.FixedAreaIntensity;
    ZProfile = conv(gaussFilter, zprofinit);
    ZProfile = ZProfile(2:end-1);
    ZProfile = ZProfile(zprofinit ~= 0);
    title_string = '';
    
    
    if isempty(zProfileFigAxes.Children)
        plot(zProfileFigAxes, currentParticleFit.z, ZProfile,'.-k');
    else
        zProfileFigAxes.Children(1).XData = currentParticleFit.z; 
        zProfileFigAxes.Children(1).YData = ZProfile;
        zProfileFigAxes.Children(1).LineStyle= '-';
        zProfileFigAxes.Children(1).Color= 'k';

    end
    
    hold(zProfileFigAxes,'on')
    
    if ~isempty(cptState.CurrentZIndex)
        if length(zProfileFigAxes.Children) < 2
            plot(zProfileFigAxes, cptState.CurrentZ, ZProfile(cptState.CurrentZIndex), 'ob')
        else
            zProfileFigAxes.Children(2).XData = cptState.CurrentZ;
            zProfileFigAxes.Children(2).YData = ZProfile(cptState.CurrentZIndex);
            zProfileFigAxes.Children(2).MarkerEdgeColor = 'b';
        end
    else
        if length(zProfileFigAxes.Children) < 2
            plot(zProfileFigAxes, cptState.CurrentZ, cptState.CurrentZ, 'or')
        else
            zProfileFigAxes.Children(2).XData = cptState.CurrentZ;
            zProfileFigAxes.Children(2).YData = cptState.CurrentZ;
            zProfileFigAxes.Children(2).MarkerEdgeColor = 'r';
        end
    end
    
    hold(zProfileFigAxes, 'off')
    
    set(zProfileFigAxes.Title, 'String', {'z-profile:';title_string});
    zProfileFigAxes.Title.FontSize = 10;
    
end

% BRIGHTEST Z-TRACE PLOT %
if ~strcmpi(ExperimentType, 'inputoutput') && ~fish
    MaxZProfile = [];
    
    CurrentParticle = cptState.CurrentParticle;
    CurrentSpots = cptState.getCurrentChannelSpots();
    CurrentParticles = cptState.getCurrentChannelParticles();
    
    % Only update the trace information if we have switched particles
    if (CurrentParticle ~= cptState.PreviousParticle) || ~exist('MaxZProfile', 'var')...
            || cptState.CurrentChannel ~= cptState.PreviousChannel
        cptState.PreviousParticle = CurrentParticle;
        PlotParticleTrace(cptState, plotTraceSettings, true);
    end
    
    for f = 1:length(cptState.Frames)
        MaxZProfile(f) = CurrentSpots(cptState.Frames(f)).Fits(...
            CurrentParticles(CurrentParticle).Index(f)).brightestZ;
    end
    
    currentFrameApproved = CurrentParticles(CurrentParticle).FrameApproved;
    if isempty(zTraceAxes.Children)
        ztrplot = plot(zTraceAxes, cptState.Frames(currentFrameApproved), MaxZProfile(currentFrameApproved), '.-k');
    else
        zTraceAxes(1).XData = cptState.Frames(currentFrameApproved);
        zTraceAxes(1).YData = MaxZProfile(currentFrameApproved);
    end
    
    hold(zTraceAxes, 'on')
    
    if length(zTraceAxes.Children) < 2
        ztrcircle = plot(zTraceAxes, cptState.Frames(cptState.Frames == cptState.CurrentFrame),...
            MaxZProfile(cptState.Frames == cptState.CurrentFrame), 'ob');
    else
        zTraceAxes(2).XData = cptState.Frames(cptState.Frames == cptState.CurrentFrame);
        zTraceAxes(2).YData = MaxZProfile(cptState.Frames == cptState.CurrentFrame);
    end
    hold(zTraceAxes, 'off')
    
    try
        xlim(zTraceAxes, [min(cptState.Frames) - 1, max(cptState.Frames) + 1]);
        ylim(zTraceAxes, [1, cptState.ZSlices + 1])
    catch
        warning('Not sure what happened here. Problem with trace fig x lim. Mention in slack if you see this, please.');
    end
else
    MaxZProfile = [];
end
end
