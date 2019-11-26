function [Frames,AmpIntegral,GaussIntegral,AmpIntegral3, ...
    ErrorIntegral, ErrorIntegral3,backGround3, ...
    AmpIntegralGauss3D, ErrorIntegralGauss3D, PreviousParticle]...
    ...
    = plotTrace(traceFigAxes, ...
    ....
    FrameInfo, CurrentChannel, PreviousChannel, ...
    CurrentParticle, PreviousParticle, lastParticle, HideApprovedFlag, lineFitted, anaphaseInMins, ...
    ElapsedTime, schnitzcells, Particles, plot3DGauss, anaphase, prophase, metaphase,prophaseInMins, metaphaseInMins,Prefix, ...
    numFrames, CurrentFrame, ZSlices, CurrentZ, Spots, ...
    correspondingNCInfo , Coefficients, ExperimentType, PreviousFrame, Frames,varargin)

%PLOTTRACE
%plot traces in checkparticletracking

if ~isempty(varargin)
    AmpIntegral = varargin{1};
    GaussIntegral = varargin{2};
    AmpIntegral3 = varargin{3};
    ErrorIntegral = varargin{4};
    ErrorIntegral3 = varargin{5};
    backGround3 = varargin{6};
    AmpIntegralGauss3D = varargin{7};
    ErrorIntegralGauss3D = varargin{8};
    FrameIndicesToFit = varargin{9};
end
if length(varargin)==10
    Spots3D = varargin{10};
end

switchParticleFlag= false;
switchFrameFlag = PreviousFrame ~= CurrentFrame;

%Only update the trace information if we have switched particles
if CurrentParticle~=PreviousParticle || ~exist('AmpIntegral', 'var')|| CurrentChannel~=PreviousChannel || lastParticle
    switchParticleFlag = true;
    switchFrameFlag = true;
    PreviousParticle=CurrentParticle;
    [Frames,AmpIntegral,GaussIntegral,AmpIntegral3,...
        ErrorIntegral, ErrorIntegral3, backGround3, ...
        AmpIntegralGauss3D, ErrorIntegralGauss3D]= ...
        PlotParticleTrace(CurrentParticle,Particles{CurrentChannel},Spots{CurrentChannel}, 'noSpline');
end


% Check if this particle has a saved manual fit or if fitInitialSlope ran
if lineFitted
    fittedXFrames = FrameIndicesToFit;
elseif  isfield(Particles{CurrentChannel},'fitApproved') && ...
        ~isempty(Particles{CurrentChannel}(CurrentParticle).fitApproved)
    approvedFit = 1;
    lineFitted = 1;
    Coefficients = Particles{CurrentChannel}(CurrentParticle).Coefficients;
    fittedXFrames = Particles{CurrentChannel}(CurrentParticle).fittedFrames;
else
    lineFitHandle = [];
    approvedFit = 0;
end

%we'll plot the spot intensity first on the left axis.
yyaxis(traceFigAxes,'left')

% finding the traceFigTimeAxis
if ~lineFitted
    traceFigTimeAxis = Frames;
else
    ncPresent = unique(correspondingNCInfo(Frames));
    priorAnaphaseInMins = anaphaseInMins(ncPresent(1)-8); %min  subtracts 8 because the first element corresponds to nc 9
    priorAnaphase = anaphase(ncPresent(1)-8); %frame
    if ~isempty(schnitzcells) && ~isempty(Particles{CurrentChannel}(CurrentParticle).Nucleus)
        nucleusFirstTimePoint = ElapsedTime(...
            schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames(1)); %min
    else
        nucleusFirstTimePoint = ElapsedTime(priorAnaphase); %min
        warning('No nucleus assigned to this particle. Using anaphase from moviedatabase as the first timepoint.')
    end
    traceFigTimeAxis = ElapsedTime(Frames) - nucleusFirstTimePoint; %min
    if exist('traceErrorBar1','var')
        delete([traceErrorBar1,traceErrorBar2,cPoint1,cPoint2])
    end
end

% plotting the lines and traces
cla(traceFigAxes)
hold(traceFigAxes, 'on')
approvedParticleFrames = Particles{CurrentChannel}(CurrentParticle).FrameApproved;
if isempty(ErrorIntegral)
    ErrorIntegral = 0;
    ErrorIntegral3 =0;
end

if plot3DGauss
    amp1 = AmpIntegral3;
    amp2 = AmpIntegralGauss3D;
    error1 = ones(length(amp1(approvedParticleFrames)),1)'*ErrorIntegral3;
    error2 = ErrorIntegralGauss3D(approvedParticleFrames);
    if lineFitted
        to = -Coefficients(2) / Coefficients(1); % minutes
        fittedXSegment = [to, traceFigTimeAxis(fittedXFrames)];
        fittedYSegment = polyval(Coefficients,fittedXSegment);
        lineFitHandle = plot(traceFigAxes,fittedXSegment,fittedYSegment);
    end
else
    amp1 = AmpIntegral;
    amp2 = AmpIntegral3;
    error1 = ones(length(amp1(approvedParticleFrames)),1)'.*ErrorIntegral';
    error2 = ones(length(amp2(approvedParticleFrames)),1)'.*ErrorIntegral3';
end

idata1 = amp1(approvedParticleFrames);
idata2 = amp2(approvedParticleFrames);

traceErrorBar1 = errorbar(traceFigAxes,traceFigTimeAxis(approvedParticleFrames),...
    idata1, error1,'.-','Color','k');
traceErrorBar2 = errorbar(traceFigAxes,traceFigTimeAxis(approvedParticleFrames),...
    idata2, error2,'.-','Color','blue');

dPoint1 = plot(traceFigAxes,traceFigTimeAxis(~approvedParticleFrames),amp1(~approvedParticleFrames),'.r');
cPoint1 = plot(traceFigAxes,traceFigTimeAxis(Frames==CurrentFrame),amp1(Frames==CurrentFrame),'ob');
dPoint2 = plot(traceFigAxes,traceFigTimeAxis(~approvedParticleFrames),amp2(~approvedParticleFrames),'.r');
cPoint2 = plot(traceFigAxes,traceFigTimeAxis(Frames==CurrentFrame),amp2(Frames==CurrentFrame),'ob');


%%
%Change the labels, xlimits, etc. only if the frame changed.
    
    % adjusting x limits
    try
        xlim(traceFigAxes,[min(traceFigTimeAxis),max(traceFigTimeAxis)]+[-1,1]);
    end
    traceFigYLimits = [0, max(AmpIntegralGauss3D)*1.1];
    
    % plotting all anaphase time points as vertical lines
    for i = 1:length(anaphase)
        if ~lineFitted
            currentAnaphaseBoundary = anaphase(i);
        else
            currentAnaphaseBoundary = anaphaseInMins(i) - priorAnaphaseInMins;
        end
        plot(traceFigAxes,ones(1,2).*currentAnaphaseBoundary,traceFigYLimits,...
            'LineWidth',2,'Color','black', 'Marker', 'none');
    end
    
    % plotting all prophase time points as vertical lines
    for i = 1:length(prophase)
        if ~lineFitted
            currentProphase = prophase(i);
        else
            currentProphase = prophaseInMins(i) - priorAnaphaseInMins;
        end
        plot(traceFigAxes,ones(1,2).*currentProphase,traceFigYLimits,...
            'LineWidth',2,'Color','blue', 'Marker', 'none');
    end
    
    % plotting all metaphase time points as vertical lines
    for i = 1:length(metaphase)
        if ~lineFitted
            currentMetaphase = metaphase(i);
        else
            currentMetaphase = metaphaseInMins(i) - priorAnaphaseInMins;
        end
        plot(traceFigAxes,ones(1,2).*currentMetaphase,traceFigYLimits,...
            'LineWidth',2,'Color','yellow');
    end
    
    % labeling plot
    ylabel(traceFigAxes,'integrated intensity (a.u.)')
    hold(traceFigAxes, 'off')
    
    % creating legend
    if plot3DGauss
        str1 = '3-slice mRNA';
        str2 = '3D-Gaussian fit mRNA';
        if ~isnan(traceFigYLimits(2))
            ylim(traceFigAxes, [0, traceFigYLimits(2)]);
        end
    else
        str1 = '1-slice mRNA';
        str2 = 'multi-slice mRNA';
    end
    
    if ~lineFitted
        legend(traceFigAxes,[traceErrorBar1,traceErrorBar2],str1,str2)
        xlabel(traceFigAxes,'frame')
    else
        legend(traceFigAxes,[traceErrorBar1,traceErrorBar2,lineFitHandle],str1,str2,...
            ['fit slope: ', num2str(round(Coefficients(1))), ' a.u./min',newline,'time on: ',num2str(roots(Coefficients)), ' min'])
        xlabel(traceFigAxes,'time since anaphase (min)')
    end
    
    if strcmpi(ExperimentType, 'inputoutput')
        yyaxis(traceFigAxes,'right')
        %now we'll plot the input protein intensity on the right-hand axis.
        plot(traceFigAxes,schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames,...
            max(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).Fluo,[],2),'r.-','DisplayName','protein')
        ylabel(traceFigAxes,'input protein intensity (a.u.)');
        hold(traceFigAxes,'on')
        plot(traceFigAxes,traceFigTimeAxis(approvedParticleFrames),AmpIntegral(~approvedParticleFrames),'.r')
        hold(traceFigAxes,'off')
    else
        traceFigAxes.YAxis(2).Visible = 'off';
    end
    
    % creating axis title
    numParticles = length(Particles{CurrentChannel});
    firstLine = [Prefix,'    Particle: ',num2str(CurrentParticle),'/',num2str(numParticles)];
    secondLine = ['Frame: ',num2str(CurrentFrame),'/',num2str(numFrames),'    ',num2str(round(FrameInfo(CurrentFrame).Time)), 's'];
    thirdLine = ['Z: ',num2str(CurrentZ),'/',num2str(ZSlices),', Ch: ',num2str(CurrentChannel)];
    
    if isfield(FrameInfo, 'nc')
        axisTitle={firstLine,...
            [secondLine,'    (nc',num2str(FrameInfo(CurrentFrame).nc),')'],...
            thirdLine};
    else
        axisTitle={firstLine,secondLine,thirdLine};
    end
    
    if HideApprovedFlag==1
        axisTitle=[axisTitle,', Showing non-flagged particles'];
    elseif HideApprovedFlag==2
        axisTitle=[axisTitle,', Showing disapproved particles'];
    end
    title(traceFigAxes,axisTitle, 'Interpreter', 'none')
    

end

