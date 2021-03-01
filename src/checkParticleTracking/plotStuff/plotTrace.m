function plotTrace(traceFigAxes, cptState, anaphaseInMins, ElapsedTime,...
    anaphase, prophase, metaphase,prophaseInMins, metaphaseInMins,Prefix,...
    numFrames, correspondingNCInfo, ExperimentType, Channels, PreProcPath, DropboxFolder,...
    plotTraceSettings)
%PLOTTRACE
%plot traces in checkparticletracking
yyaxis(traceFigAxes,'left');
if cptState.UseTwinTraces & ~isempty(cptState.TwinParticle)
    switchTwinParticleFlag= false;
    
    % Only update the trace information if we have switched particles
    
    switchParticleFlag = true;
    switchFrameFlag = true;
    cptState.PreviousTwinParticle = cptState.TwinParticle;
    
    PlotParticleTrace(cptState, plotTraceSettings, true, true);
    
    % Check if this particle has a saved manual fit or if fitInitialSlope ran
    if cptState.lineFitted
        fittedXFrames = cptState.FrameIndicesToFit;
    elseif  isfield(cptState.Particles{cptState.CurrentChannelIndex},'fitApproved') && ...
            ~isempty(cptState.getTwinParticle().fitApproved)
        approvedFit = 1;
        
        % JP: are this changes to lineFitted and Coefficients supposed to be visible
        % from outside this function? (ie. change it on global cptState property?)
        cptState.lineFitted = 1;
        cptState.Coefficients = cptState.getTwinParticle().Coefficients;
        fittedXFrames = cptState.getTwinParticle().fittedFrames;
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
        if ~isempty(cptState.schnitzcells) && ~isempty(cptState.getTwinParticle().Nucleus)
            nucleusFirstTimePoint = ElapsedTime(...
                cptState.schnitzcells(cptState.getTwinParticle().Nucleus).frames(1)); %min
        else
            nucleusFirstTimePoint = ElapsedTime(priorAnaphase); %min
            warning('No nucleus assigned to this particle. Using anaphase from moviedatabase as the first timepoint.')
        end
        traceFigTimeAxis = ElapsedTime(cptState.Frames) - nucleusFirstTimePoint; %min
        if exist('cPoint3','var')
            delete([traceErrorBar3, traceErrorBar4, cPoint3,cPoint4])
        end
    end
    
    % plotting the lines and traces
    hold(traceFigAxes, 'on')
    approvedParticleFrames = cptState.getTwinParticle().FrameApproved;
    if isempty(plotTraceSettings.ErrorIntegral)
        plotTraceSettings.ErrorIntegral = 0;
        plotTraceSettings.ErrorIntegral3 = 0;
    end
    
    if cptState.plot3DGauss
        amp1 = plotTraceSettings.AmpIntegral3;
        amp2 = plotTraceSettings.AmpIntegralGauss3D;
        
        error1aux = plotTraceSettings.ErrorIntegral3;
        error1 = ones(length(amp1(approvedParticleFrames)),1)'*error1aux;
        
        error2 = plotTraceSettings.ErrorIntegralGauss3D(approvedParticleFrames);
        
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
    traceErrorBar3 = traceFigAxes.Children(end-2);
    traceErrorBar4 = traceFigAxes.Children(end-3);
    set(traceErrorBar1, 'HandleVisibility', 'off');
    set(traceErrorBar2, 'HandleVisibility', 'off');
    
    
    set(traceErrorBar3, 'XData', traceFigTimeAxis(approvedParticleFrames),...
        'YData', idata1,'YNegativeDelta', error1, 'YPositiveDelta', error1, 'HandleVisibility', 'off');
    set(traceErrorBar4, 'XData', traceFigTimeAxis(approvedParticleFrames),...
        'YData', idata2,'YNegativeDelta', error2, 'YPositiveDelta', error2, 'HandleVisibility', 'off');
    
    cla(traceFigAxes);
    
    
    set(traceErrorBar4, 'HandleVisibility', 'on', 'Visible', 'on');
    set(traceErrorBar3, 'HandleVisibility', 'on', 'Visible', 'on');
    
    dPoint3 = plot(traceFigAxes,traceFigTimeAxis(~approvedParticleFrames),amp1(~approvedParticleFrames),'.r');
    cPoint3 = plot(traceFigAxes,traceFigTimeAxis(cptState.Frames==cptState.CurrentFrame),amp1(cptState.Frames==cptState.CurrentFrame),'og');
    dPoint4 = plot(traceFigAxes,traceFigTimeAxis(~approvedParticleFrames),amp2(~approvedParticleFrames),'.r');
    cPoint4 = plot(traceFigAxes,traceFigTimeAxis(cptState.Frames==cptState.CurrentFrame),amp2(cptState.Frames==cptState.CurrentFrame),'om');
    
    
    set(dPoint3, 'HandleVisibility', 'on', 'Visible', 'on');
    set(cPoint3, 'HandleVisibility', 'on', 'Visible', 'on');
    set(dPoint4, 'HandleVisibility', 'on', 'Visible', 'on');
    set(cPoint4, 'HandleVisibility', 'on', 'Visible', 'on');
    
elseif cptState.UseTwinTraces & isempty(cptState.TwinParticle)
    traceErrorBar1 = traceFigAxes.Children(end);
    traceErrorBar2 = traceFigAxes.Children(end-1);
    traceErrorBar3 = traceFigAxes.Children(end-2);
    traceErrorBar4 = traceFigAxes.Children(end-3);
    set(traceErrorBar1, 'HandleVisibility', 'off');
    set(traceErrorBar2, 'HandleVisibility', 'off');
    
    
    set(traceErrorBar3, 'HandleVisibility', 'off');
    set(traceErrorBar4, 'HandleVisibility', 'off');
    
    cla(traceFigAxes);
    
    set(traceErrorBar3, 'Visible', 'off');
    set(traceErrorBar4, 'Visible', 'off');
    
    
end

%%
switchParticleFlag= false;
switchFrameFlag = cptState.PreviousFrame ~= cptState.CurrentFrame;

% Only update the trace information if we have switched particles
if cptState.UseTwinTraces || cptState.CurrentParticle ~= cptState.PreviousParticle || isempty(plotTraceSettings.AmpIntegral) || cptState.CurrentChannelIndex ~= cptState.PreviousChannel || cptState.lastParticle
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
    if ~isempty(cptState.schnitzcells) && ~isempty(cptState.getCurrentParticle().Nucleus)
        nucleusFirstTimePoint = ElapsedTime(...
            cptState.schnitzcells(cptState.getCurrentParticle().Nucleus).frames(1)); %min
    else
        nucleusFirstTimePoint = ElapsedTime(priorAnaphase); %min
        warning('No nucleus assigned to this particle. Using anaphase from moviedatabase as the first timepoint.')
    end
    traceFigTimeAxis = ElapsedTime(cptState.Frames) - nucleusFirstTimePoint; %min
    if exist('cPoint1','var')
        delete([cPoint1,cPoint2])
    end
end

% plotting the lines and traces
hold(traceFigAxes, 'on')
approvedParticleFrames = cptState.getCurrentParticle().FrameApproved;
if isempty(plotTraceSettings.ErrorIntegral)
    plotTraceSettings.ErrorIntegral = 0;
    plotTraceSettings.ErrorIntegral3 = 0;
end

if cptState.plot3DGauss
    amp1 = plotTraceSettings.AmpIntegral3;
    amp2 = plotTraceSettings.AmpIntegralGauss3D;
    
    error1aux = plotTraceSettings.ErrorIntegral3;
    error1 = ones(length(amp1(approvedParticleFrames)),1)'*error1aux;
    
    error2 = plotTraceSettings.ErrorIntegralGauss3D(approvedParticleFrames);
    
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



if ~cptState.UseTwinTraces
    traceErrorBar1 = traceFigAxes.Children(end);
    traceErrorBar2 = traceFigAxes.Children(end-1);
end
set(traceErrorBar1, 'XData', traceFigTimeAxis(approvedParticleFrames),...
    'YData', idata1,'YNegativeDelta', error1, 'YPositiveDelta', error1, 'HandleVisibility', 'off');
set(traceErrorBar2, 'XData', traceFigTimeAxis(approvedParticleFrames),...
    'YData', idata2,'YNegativeDelta', error2, 'YPositiveDelta', error2, 'HandleVisibility', 'off');


if ~cptState.UseTwinTraces
    traceErrorBar3 = traceFigAxes.Children(end-2);
    traceErrorBar4 = traceFigAxes.Children(end-3);
    set(traceErrorBar4, 'HandleVisibility', 'off');
    set(traceErrorBar3, 'HandleVisibility', 'off');
    
    cla(traceFigAxes);
    set(traceErrorBar4, 'HandleVisibility', 'on');
    set(traceErrorBar3, 'HandleVisibility', 'on');
end
set(traceErrorBar2, 'HandleVisibility', 'on', 'Visible', 'on');
set(traceErrorBar1, 'HandleVisibility', 'on', 'Visible', 'on');

dPoint1 = plot(traceFigAxes,traceFigTimeAxis(~approvedParticleFrames),amp1(~approvedParticleFrames),'.r');
cPoint1 = plot(traceFigAxes,traceFigTimeAxis(cptState.Frames==cptState.CurrentFrame),amp1(cptState.Frames==cptState.CurrentFrame),'ob');
dPoint2 = plot(traceFigAxes,traceFigTimeAxis(~approvedParticleFrames),amp2(~approvedParticleFrames),'.r');
cPoint2 = plot(traceFigAxes,traceFigTimeAxis(cptState.Frames==cptState.CurrentFrame),amp2(cptState.Frames==cptState.CurrentFrame),'ob');

set(dPoint1, 'HandleVisibility', 'on', 'Visible', 'on');
set(cPoint1, 'HandleVisibility', 'on', 'Visible', 'on');
set(dPoint2, 'HandleVisibility', 'on', 'Visible', 'on');
set(cPoint2, 'HandleVisibility', 'on', 'Visible', 'on');
if isempty(traceFigTimeAxis(~approvedParticleFrames))
    if exist('dPoint3', 'var')
        set(dPoint3, 'Visible', 'off');
    end
    if exist('dPoint4', 'var')
        set(dPoint4, 'Visible', 'off');
    end
    if exist('cPoint3', 'var')
        set(cPoint3, 'Visible', 'off');
    end
    if exist('cPoint4', 'var')
        set(cPoint4, 'Visible', 'off');
    end
end
%%
%Change the labels, xlimits, etc. only if the frame changed.
ncPresent = unique(correspondingNCInfo(cptState.Frames));
priorAnaphaseInMins = anaphaseInMins(ncPresent(1)-8); %min  subtracts 8 because the first element corresponds to nc 9
priorAnaphase = anaphase(ncPresent(1)-8); %frame
if ncPresent(1) == 14
    subsequentAnaphase = numFrames;
else
    subsequentAnaphase =  anaphase(ncPresent(1)-7)-1;
end
DataXmaxes = NaN(1, 4);
DataXmins= NaN(1, 4);
DataYmaxes = NaN(1,4);
DataYmins = NaN(1,4);
try
    XData1 = get(traceErrorBar1, 'XData');
    DataXmins(1) = nanmin(XData1);
    DataXmaxes(1) =nanmax(XData1);
    YData1 = get(traceErrorBar1, 'YData');
    DataYmins(1) = nanmin(YData1);
    DataYmaxes(1) =nanmax(YData1);
end
try
    XData2 = get(traceErrorBar2, 'XData');
    DataXmins(2) = nanmin(XData2);
    DataXmaxes(2) =nanmax(XData2);
    YData2 = get(traceErrorBar2, 'YData');
    DataYmins(2) = nanmin(YData2);
    DataYmaxes(2) =nanmax(YData2);
end
if cptState.UseTwinTraces & ~isempty(cptState.TwinParticle)
    try
        XData3 = get(traceErrorBar3, 'XData');
        DataXmins(3) = nanmin(XData3);
        DataXmaxes(3) =nanmax(XData3);
        YData3 = get(traceErrorBar3, 'YData');
        DataYmins(3) = nanmin(YData3);
        DataYmaxes(3) =nanmax(YData3);
    end
    try
        XData4 = get(traceErrorBar4, 'XData');
        DataXmins(4) = nanmin(XData4);
        DataXmaxes(4) =nanmax(XData4);
        YData4 = get(traceErrorBar4, 'YData');
        DataYmins(4) = nanmin(YData4);
        DataYmaxes(4) =nanmax(YData4);
    end
end
XlimMin = min(nanmin(DataXmins), priorAnaphase);
XlimMax = max(nanmax(DataXmaxes), subsequentAnaphase);
XlimMax = nanmax(DataXmaxes);
YlimMax = nanmax(DataYmaxes);
try
    setPlotsInvisible(traceFigAxes);
    xlim(traceFigAxes,[XlimMin, XlimMax]+[-1,1]);
    setPlotsVisible(traceFigAxes);
end
traceFigYLimits = [0, nanmax(YlimMax, 1)*1.3];
ylim(traceFigAxes,traceFigYLimits);

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

if strcmpi(ExperimentType, 'inputoutput') | (cptState.PlotInputChannel & ~isempty(cptState.liveExperiment.inputChannels))
    yyaxis(traceFigAxes,'right')
    %now we'll plot the input protein intensity on the right-hand axis.
    if ~isfield(cptState.schnitzcells, 'Fluo')
        dummy = cell(length(cptState.schnitzcells), 1);
        [cptState.schnitzcells.Fluo] = dummy{:};
    else
        proteinLine = traceFigAxes.Children(end);
        if ~isempty(cptState.getCurrentParticle().Nucleus)
            proteinFluo = cptState.schnitzcells(cptState.getCurrentParticle().Nucleus).Fluo;
            if ~isempty(proteinFluo)
                set(proteinLine, 'XData', cptState.schnitzcells(cptState.getCurrentParticle().Nucleus).frames,...
                    'YData', max(proteinFluo,[],2),  'Visible', 'on');
            else
                warning('protein fluo empty. maybe rerun integrateschnitzfluo?');
            end
            try
                setPlotsInvisible(traceFigAxes);
                pXdata = get(proteinLine, 'XData');
                xlim(traceFigAxes,double([nanmin(XlimMin, nanmin(pXdata)), nanmax(XlimMax, nanmax(pXdata))])+[-1,1]);
                setPlotsVisible(traceFigAxes);
            end
        else
            set(proteinLine, 'Visible', 'off');
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


if cptState.HideSingleSliceTrace
    set(traceErrorBar1, 'Visible', 'off');
    set(get(get(traceErrorBar1, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    set(cPoint1, 'Visible', 'off');
    set(dPoint1, 'Visible', 'off');
    if cptState.UseTwinTraces & ~isempty(cptState.TwinParticle)
        set(traceErrorBar3, 'Visible', 'off');
        set(get(get(traceErrorBar3, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        set(cPoint3, 'Visible', 'off');
        set(dPoint3, 'Visible', 'off');
    end
else
    set(traceErrorBar1, 'Visible', 'on');
    set(get(get(traceErrorBar1, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
    set(cPoint1, 'Visible', 'on');
    set(dPoint1, 'Visible', 'on');
    if cptState.UseTwinTraces
        set(traceErrorBar3, 'Visible', 'on');
        set(get(get(traceErrorBar3, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'on');
        set(cPoint3, 'Visible', 'on');
        set(dPoint3, 'Visible', 'on');
    end
end

if cptState.UseTwinTraces
    set(traceErrorBar3, 'HandleVisibility', 'on');
    set(traceErrorBar4, 'HandleVisibility', 'on');
end


end