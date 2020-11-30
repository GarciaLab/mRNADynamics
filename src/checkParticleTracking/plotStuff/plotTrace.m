function plotTrace(traceFigAxes, cptState, anaphaseInMins, ElapsedTime,...
    anaphase, prophase, metaphase,prophaseInMins, metaphaseInMins,Prefix,...
    numFrames, correspondingNCInfo, ExperimentType, Channels, PreProcPath, DropboxFolder,...
    plotTraceSettings)
%PLOTTRACE
%plot traces in checkparticletracking

switchParticleFlag= false;
switchFrameFlag = cptState.PreviousFrame ~= cptState.CurrentFrame;

% Only update the trace information if we have switched particles
if cptState.CurrentParticle ~= cptState.PreviousParticle || isempty(plotTraceSettings.AmpIntegral) || cptState.CurrentChannelIndex ~= cptState.PreviousChannel || cptState.lastParticle
    switchParticleFlag = true;
    switchFrameFlag = true;
    cptState.PreviousParticle = cptState.CurrentParticle;
    PlotParticleTrace(cptState, plotTraceSettings, true);
end


% Check if this particle has a saved manual fit or if fitInitialSlope ran
if cptState.lineFitted
    fittedXFrames = cptState.FrameIndicesToFit;
elseif  isfield(cptState.Particles{cptState.CurrentChannelIndex},'fitApproved') && ...
        ~isempty(cptState.getCurrentParticle().fitApproved)
    approvedFit = 1;

    % JP: are this changes to lineFitted and Coefficients supposed to be visible
    % from outside this function? (ie. change it on global cptState property?)
    cptState.lineFitted = 1;
    cptState.Coefficients = cptState.getCurrentParticle().Coefficients;
    fittedXFrames = cptState.getCurrentParticle().fittedFrames;
else
    lineFitHandle = [];
    approvedFit = 0;
end

% we'll plot the spot intensity first on the left axis.
yyaxis(traceFigAxes,'left')

% finding the traceFigTimeAxis
if ~cptState.lineFitted
    traceFigTimeAxis = cptState.Frames;
else
    ncPresent = unique(correspondingNCInfo(cptState.Frames));
    priorAnaphaseInMins = anaphaseInMins(ncPresent(1)-8); %min  subtracts 8 because the first element corresponds to nc 9
    priorAnaphase = anaphase(ncPresent(1)-8); %frame
    if ~isempty(cptState.schnitzcells) && ~isempty(cptStategetCurrentParticle().Nucleus)
        nucleusFirstTimePoint = ElapsedTime(...
            cptState.schnitzcells(cptState.getCurrentParticle().Nucleus).frames(1)); %min
    else
        nucleusFirstTimePoint = ElapsedTime(priorAnaphase); %min
        warning('No nucleus assigned to this particle. Using anaphase from moviedatabase as the first timepoint.')
    end
    traceFigTimeAxis = ElapsedTime(cptState.Frames) - nucleusFirstTimePoint; %min
    if exist('traceErrorBar1','var')
        delete([traceErrorBar1,traceErrorBar2,cPoint1,cPoint2])
    end
end

% plotting the lines and traces
hold(traceFigAxes, 'on')
approvedParticleFrames = cptState.getCurrentParticle().FrameApproved;
if isempty(plotTraceSettings.ErrorIntegral)
    plotTraceSettings.ErrorIntegral = 0;
    plotTraceSettings.ErrorIntegral3 = 0;
end

lgd_str = {'mRNA (2D 1 slice)'};
if cptState.plot3DGauss
    amp1 = plotTraceSettings.AmpIntegral;
    amp2 = plotTraceSettings.gauss3DIntensity;        
    error1aux = plotTraceSettings.ErrorIntegral3;
    error1 = ones(length(amp1(approvedParticleFrames)),1)'*error1aux;
    
    error2 = plotTraceSettings.gauss3DIntensitySE(approvedParticleFrames);
    
    if cptState.lineFitted
        to = -cptState.Coefficients(2) / cptState.Coefficients(1); % minutes
        fittedXSegment = [to, traceFigTimeAxis(fittedXFrames)];
        fittedYSegment = polyval(cptState.Coefficients,fittedXSegment);
        lineFitHandle = plot(traceFigAxes,fittedXSegment,fittedYSegment);
    end
else
    amp1 = plotTraceSettings.AmpIntegral;
    amp2 = plotTraceSettings.AmpIntegral3;    
    error1aux = plotTraceSettings.ErrorIntegral;
    error1 = ones(length(amp1(approvedParticleFrames)),1)'.*error1aux';

    error2aux = plotTraceSettings.ErrorIntegral3;
    error2 = ones(length(amp2(approvedParticleFrames)),1)'.*error2aux';
end

idata1 = amp1(approvedParticleFrames);
idata2 = amp2(approvedParticleFrames);


traceErrorBar1 = traceFigAxes.Children(end);
traceErrorBar2 = traceFigAxes.Children(end-1);
set(traceErrorBar1, 'XData', traceFigTimeAxis(approvedParticleFrames),...
    'YData', idata1,'YNegativeDelta', error1, 'YPositiveDelta', error1, 'HandleVisibility', 'off');
set(traceErrorBar2, 'XData', traceFigTimeAxis(approvedParticleFrames),...
    'YData', idata2,'YNegativeDelta', error2, 'YPositiveDelta', error2, 'HandleVisibility', 'off');

cla(traceFigAxes);
set(traceErrorBar2, 'HandleVisibility', 'on');
set(traceErrorBar1, 'HandleVisibility', 'on');

dPoint1 = plot(traceFigAxes,traceFigTimeAxis(~approvedParticleFrames),amp1(~approvedParticleFrames),'.r');
cPoint1 = plot(traceFigAxes,traceFigTimeAxis(cptState.Frames==cptState.CurrentFrame),amp1(cptState.Frames==cptState.CurrentFrame),'ob');
dPoint2 = plot(traceFigAxes,traceFigTimeAxis(~approvedParticleFrames),amp2(~approvedParticleFrames),'.r');
cPoint2 = plot(traceFigAxes,traceFigTimeAxis(cptState.Frames==cptState.CurrentFrame),amp2(cptState.Frames==cptState.CurrentFrame),'ob');


%%
%Change the labels, xlimits, etc. only if the frame changed.
    
    % adjusting x limits
    try
        setPlotsInvisible(traceFigAxes);
        xlim(traceFigAxes,[min(traceFigTimeAxis),max(traceFigTimeAxis)]+[-1,1]);
        setPlotsVisible(traceFigAxes);
    end
    if cptState.plot3DGauss
        traceFigYLimits = [0, max(plotTraceSettings.gauss3DIntensity)*1.1];
    else
        traceFigYLimits = [0, max(plotTraceSettings.AmpIntegral3)*1.1];
    end
    
    % plotting all anaphase time points as vertical lines
    for i = 1:length(anaphase)
        if ~cptState.lineFitted
            currentAnaphaseBoundary = anaphase(i);
        else
            currentAnaphaseBoundary = anaphaseInMins(i) - priorAnaphaseInMins;
        end
        plot(traceFigAxes,ones(1,2).*double(currentAnaphaseBoundary),traceFigYLimits,...
            'LineWidth',2,'Color','black', 'Marker', 'none');
    end
    
    % plotting all prophase time points as vertical lines
    for i = 1:length(prophase)
        if ~cptState.lineFitted
            currentProphase = prophase(i);
        else
            currentProphase = prophaseInMins(i) - priorAnaphaseInMins;
        end
        plot(traceFigAxes,ones(1,2).*currentProphase,traceFigYLimits,...
            'LineWidth',2,'Color','blue', 'Marker', 'none');
    end
    
    % plotting all metaphase time points as vertical lines
    for i = 1:length(metaphase)
        if ~cptState.lineFitted
            currentMetaphase = metaphase(i);
        else
            currentMetaphase = metaphaseInMins(i) - priorAnaphaseInMins;
        end
        plot(traceFigAxes,ones(1,2).*currentMetaphase,traceFigYLimits,...
            'LineWidth',2,'Color','yellow');
    end
    
    hold(traceFigAxes, 'off')
    
    % creating legend
    if cptState.plot3DGauss && ~isnan(traceFigYLimits(2))
            setPlotsInvisible(traceFigAxes);
            ylim(traceFigAxes, [0, traceFigYLimits(2)]);
            setPlotsVisible(traceFigAxes);
    end

    
    if cptState.lineFitted
        legend(traceFigAxes,[traceErrorBar1,traceErrorBar2,lineFitHandle],str1,str2,...
            ['fit slope: ', num2str(round(cptState.Coefficients(1))), ' a.u./min',newline,'time on: ',num2str(roots(cptState.Coefficients)), ' min'])
        xlabel(traceFigAxes,'time since anaphase (min)')
    end
    
      hold(traceFigAxes,'on')
      plot(traceFigAxes,traceFigTimeAxis(~approvedParticleFrames),plotTraceSettings.AmpIntegral(~approvedParticleFrames),'.r')
      hold(traceFigAxes,'off')
    
    if strcmpi(ExperimentType, 'inputoutput')        
        yyaxis(traceFigAxes,'right')
        %now we'll plot the input protein intensity on the right-hand axis.
        if ~isfield(cptState.schnitzcells, 'Fluo')
                dummy = cell(length(cptState.schnitzcells), 1);
                [cptState.schnitzcells.Fluo] = dummy{:};
        else
            proteinLine = traceFigAxes.Children(end);
            proteinFluo = cptState.schnitzcells(cptState.getCurrentParticle().Nucleus).Fluo;
            if ~isempty(proteinFluo)
                set(proteinLine, 'XData', cptState.schnitzcells(cptState.getCurrentParticle().Nucleus).frames,...
                    'YData', max(proteinFluo,[],2));                
            else
                warning('protein fluo empty. maybe rerun integrateschnitzfluo?');
            end
        end
    end 
    
    % creating axis title
    numParticles = length(cptState.Particles{cptState.CurrentChannelIndex});
    firstLine = [Prefix,'    Particle: ',num2str(cptState.CurrentParticle),'/',num2str(numParticles)];
    secondLine = ['Frame: ',num2str(cptState.CurrentFrame),'/',num2str(numFrames),'    ',num2str(round(cptState.FrameInfo(cptState.CurrentFrame).Time)), 's'];
    thirdLine = ['Z: ',num2str(cptState.CurrentZ),'/',num2str(cptState.ZSlices),', Ch: ',num2str(cptState.CurrentChannelIndex)];
    
    if isfield(cptState.FrameInfo, 'nc')
        axisTitle={firstLine,...
            [secondLine,'    (nc',num2str(cptState.FrameInfo(cptState.CurrentFrame).nc),')'],...
            thirdLine};
    else
        axisTitle={firstLine,secondLine,thirdLine};
    end
    
    if cptState.HideApprovedFlag == 1
        axisTitle=[axisTitle,', Showing non-flagged particles'];
    elseif cptState.HideApprovedFlag == 2
        axisTitle=[axisTitle,', Showing disapproved particles'];
    end
    hold(traceFigAxes, 'off');
    setPlotsInvisible(traceFigAxes);
    set(traceFigAxes.Title,'String', axisTitle);
    setPlotsVisible(traceFigAxes);

end

