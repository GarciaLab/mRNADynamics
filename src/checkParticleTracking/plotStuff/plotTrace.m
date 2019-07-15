function [Frames,AmpIntegral,GaussIntegral,AmpIntegral3,AmpIntegral5, ...
    ErrorIntegral, ErrorIntegral3, ErrorIntegral5,backGround3, ...
    AmpIntegralGauss3D, ErrorIntegralGauss3D, PreviousParticle]...
    ...
    = plotTrace(traceFigAxes, ...
    ....
    FrameInfo, CurrentChannel, PreviousChannel, ...
    CurrentParticle, PreviousParticle, lastParticle, HideApprovedFlag, lineFitted, anaphaseInMins, ...
    ElapsedTime, schnitzcells, Particles, plot3DGauss, anaphase, prophase, metaphase,prophaseInMins, metaphaseInMins,Prefix, ...
    numFrames, CurrentFrame, ZSlices, CurrentZ, Spots, ...
    correspondingNCInfo , Coefficients, ExperimentType, ...
    Frames, AmpIntegral, GaussIntegral, AmpIntegral3, AmpIntegral5, ...
    ErrorIntegral, ErrorIntegral3, ErrorIntegral5, backGround3, ...
    AmpIntegralGauss3D, ErrorIntegralGauss3D,FrameIndicesToFit)
%PLOTTRACE Summary of this function goes here
%   Detailed explanation goes here

switchParticleFlag= false;
switchFrameFlag = false;

%Only update the trace information if we have switched particles
if (CurrentParticle~=PreviousParticle)||~exist('AmpIntegral', 'var')||(CurrentChannel~=PreviousChannel) || lastParticle
    switchParticleFlag = true;
    PreviousParticle=CurrentParticle;
    [Frames,AmpIntegral,GaussIntegral,AmpIntegral3,AmpIntegral5, ...
        ErrorIntegral, ErrorIntegral3, ErrorIntegral5,backGround3, ...
        AmpIntegralGauss3D, ErrorIntegralGauss3D]= ...
        PlotParticleTrace(CurrentParticle,Particles{CurrentChannel},Spots{CurrentChannel});
end
if switchParticleFlag
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
    if ~plot3DGauss
        traceErrorBar1 = errorbar(traceFigAxes, traceFigTimeAxis(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
            AmpIntegral(Particles{CurrentChannel}(CurrentParticle).FrameApproved),ones(length(AmpIntegral(Particles{CurrentChannel}(CurrentParticle).FrameApproved)),1)'*ErrorIntegral,'.-k');
        traceErrorBar2 = errorbar(traceFigAxes,traceFigTimeAxis(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
            AmpIntegral3(Particles{CurrentChannel}(CurrentParticle).FrameApproved),ones(length(AmpIntegral3(Particles{CurrentChannel}(CurrentParticle).FrameApproved)),1)'*ErrorIntegral3,'.-','Color','green');
        dPoint1 = plot(traceFigAxes,traceFigTimeAxis(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),AmpIntegral(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.r');
        cPoint1 = plot(traceFigAxes,traceFigTimeAxis(Frames==CurrentFrame),AmpIntegral(Frames==CurrentFrame),'ob');
        dPoint2 = plot(traceFigAxes,traceFigTimeAxis(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),AmpIntegral3(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.r');
        cPoint2 = plot(traceFigAxes,traceFigTimeAxis(Frames==CurrentFrame),AmpIntegral3(Frames==CurrentFrame),'ob');
    elseif lineFitted
        % plotting the traces
        traceErrorBar1 = errorbar(traceFigAxes,traceFigTimeAxis(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
            AmpIntegral3(Particles{CurrentChannel}(CurrentParticle).FrameApproved),ones(length(AmpIntegral3(Particles{CurrentChannel}(CurrentParticle).FrameApproved)),1)'*ErrorIntegral3,'.-','Color','green');
        traceErrorBar2 = plot(traceFigAxes,traceFigTimeAxis(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
            AmpIntegralGauss3D(Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.-','Color','blue');

        % calculate the fittedXSegment and fittedYSegment
        to = -Coefficients(2) / Coefficients(1); % minutes 
        fittedXSegment = [to, traceFigTimeAxis(fittedXFrames)];
        fittedYSegment = polyval(Coefficients,fittedXSegment);
        lineFitHandle = plot(traceFigAxes,fittedXSegment,fittedYSegment); 

        dPoint1 = plot(traceFigAxes,traceFigTimeAxis(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),AmpIntegral(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.r');
        cPoint1 = plot(traceFigAxes,traceFigTimeAxis(Frames==CurrentFrame),AmpIntegral3(Frames==CurrentFrame),'ob');
        dPoint2 = plot(traceFigAxes,traceFigTimeAxis(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),AmpIntegral3(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.r');
        cPoint2 = plot(traceFigAxes,traceFigTimeAxis(Frames==CurrentFrame),AmpIntegralGauss3D(Frames==CurrentFrame),'ob');
    else
        traceErrorBar1 = errorbar(traceFigAxes,traceFigTimeAxis(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
            AmpIntegral3(Particles{CurrentChannel}(CurrentParticle).FrameApproved),ones(length(AmpIntegral3(Particles{CurrentChannel}(CurrentParticle).FrameApproved)),1)'*ErrorIntegral3,'.-','Color','green');
    %     traceErrorBar2 = plot(traceFigAxes,traceFigTimeAxis(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
    %         AmpIntegralGauss3D(Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.-','Color','blue');
        traceErrorBar2 = errorbar(traceFigAxes,traceFigTimeAxis(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
            AmpIntegralGauss3D(Particles{CurrentChannel}(CurrentParticle).FrameApproved),ErrorIntegralGauss3D(Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.-','Color','blue');     

        dPoint1 = plot(traceFigAxes,traceFigTimeAxis(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),AmpIntegral(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.r');
        cPoint1 = plot(traceFigAxes,traceFigTimeAxis(Frames==CurrentFrame),AmpIntegral3(Frames==CurrentFrame),'ob');
        dPoint2 = plot(traceFigAxes,traceFigTimeAxis(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),AmpIntegral3(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.r');
        cPoint2 = plot(traceFigAxes,traceFigTimeAxis(Frames==CurrentFrame),AmpIntegralGauss3D(Frames==CurrentFrame),'ob');
    end

    % adjusting x limits
    try
        xlim(traceFigAxes,[min(traceFigTimeAxis),max(traceFigTimeAxis)]+[-1,1]);
    catch
        %             error('Not sure what happened here. Problem with trace fig x lim. Talk to AR if you see this, please.');
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
        plot(traceFigAxes,Frames(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),AmpIntegral(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.r')
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

end

