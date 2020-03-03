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
        
        MaxZ = currentParticleFit.brightestZ;
        
        plot(zProfileFigAxes, currentParticleFit.z, ZProfile,'.-k');
        
        hold(zProfileFigAxes,'on')
        if ~isempty(cptState.CurrentZIndex)
            plot(zProfileFigAxes, cptState.CurrentZ, ZProfile(cptState.CurrentZIndex), 'ob')
        else
            plot(zProfileFigAxes, cptState.CurrentZ, cptState.CurrentZ, 'or')
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
        if (CurrentParticle ~= cptState.PreviousParticle) || ~exist('MaxZProfile', 'var') || cptState.CurrentChannel ~= cptState.PreviousChannel
            cptState.PreviousParticle = CurrentParticle;
            PlotParticleTrace(cptState, plotTraceSettings, true);
        end

        for i = 1:length(cptState.Frames)
            MaxZProfile(i) = CurrentSpots(cptState.Frames(i)).Fits(CurrentParticles(CurrentParticle).Index(i)).brightestZ;
        end

        currentFrameApproved = CurrentParticles(CurrentParticle).FrameApproved;
        plot(zTraceAxes, cptState.Frames(currentFrameApproved), MaxZProfile(currentFrameApproved), '.-k');
        hold(zTraceAxes, 'on')
        plot(zTraceAxes, cptState.Frames(cptState.Frames == cptState.CurrentFrame), MaxZProfile(cptState.Frames == cptState.CurrentFrame), 'ob');
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
