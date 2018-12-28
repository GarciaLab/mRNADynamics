function [Frames,AmpIntegral,GaussIntegral,AmpIntegral3,AmpIntegral5, ...
            ErrorIntegral, ErrorIntegral3, ErrorIntegral5,backGround3, ...
            AmpIntegralGauss3D, ErrorIntegralGauss3D, PreviousParticle] = plotTrace(traceFigAxes, ...
            FrameInfo, CurrentChannel, PreviousChannel, ...
    CurrentParticle, PreviousParticle, lastParticle, HideApprovedFlag, lineFit, anaphaseInMins,...
    ElapsedTime, schnitzcells, Particles, plot3DGauss, anaphase,prophase, metaphase, prophaseInMins, metaphaseInMins,Prefix, ...
    DefaultDropboxFolder, numFrames, CurrentFrame, ZSlices, CurrentZ, Spots, ...
    correspondingNCInfo, fit1E, Coefficients, Frames, AmpIntegral, GaussIntegral, AmpIntegral3, ...
    AmpIntegral5, ErrorIntegral, ErrorIntegral3, ErrorIntegral5, backGround3, ...
    AmpIntegralGauss3D, ErrorIntegralGauss3D)
%PLOTTRACE Summary of this function goes here
%   Detailed explanation goes here

numParticles = length(Particles{CurrentChannel});

%Only update the trace information if we have switched particles
if (CurrentParticle~=PreviousParticle)||~exist('AmpIntegral', 'var')||(CurrentChannel~=PreviousChannel) || lastParticle
    PreviousParticle=CurrentParticle;
    [Frames,AmpIntegral,GaussIntegral,AmpIntegral3,AmpIntegral5, ...
        ErrorIntegral, ErrorIntegral3, ErrorIntegral5,backGround3, ...
        AmpIntegralGauss3D, ErrorIntegralGauss3D]= ...
        PlotParticleTrace(CurrentParticle,Particles{CurrentChannel},Spots{CurrentChannel});
end
%         yyaxis(traceFigAxes,'left');
if ~lineFit
    traceFigTimeAxis = Frames;
    cla(traceFigAxes)
else
    ncPresent = unique(correspondingNCInfo(Frames));
    % below subtracts 8 because the first element corresponds to nc 9
    priorAnaphaseInMins = anaphaseInMins(ncPresent(1)-8);
    nucleusFirstFrame = ElapsedTime(...
        schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames(1));
%             traceFigTimeAxis = ElapsedTime(Frames);
%             traceFigTimeAxis = ElapsedTime(Frames) - priorAnaphaseInMins;
    traceFigTimeAxis = ElapsedTime(Frames) - nucleusFirstFrame;
    if exist('traceErrorBar1','var')
        delete([traceErrorBar1,traceErrorBar2,cPoint1,cPoint2])
    end
end

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

else
   traceErrorBar1 = errorbar(traceFigAxes,traceFigTimeAxis(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
        AmpIntegral3(Particles{CurrentChannel}(CurrentParticle).FrameApproved),ones(length(AmpIntegral3(Particles{CurrentChannel}(CurrentParticle).FrameApproved)),1)'*ErrorIntegral3,'.-','Color','green');
   traceErrorBar2 = plot(traceFigAxes,Frames(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
            AmpIntegralGauss3D(Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.-','Color','blue');
    dPoint1 = plot(traceFigAxes,traceFigTimeAxis(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),AmpIntegral(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.r');
    cPoint1 = plot(traceFigAxes,traceFigTimeAxis(Frames==CurrentFrame),AmpIntegral3(Frames==CurrentFrame),'ob');
    dPoint2 = plot(traceFigAxes,traceFigTimeAxis(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),AmpIntegral3(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.r');
    cPoint2 = plot(traceFigAxes,traceFigTimeAxis(Frames==CurrentFrame),AmpIntegralGauss3D(Frames==CurrentFrame),'ob');

end
    %         p3 = errorbar(traceFigAxes,Frames(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
%            AmpIntegral5(Particles{CurrentChannel}(CurrentParticle).FrameApproved),ones(length(AmpIntegral5(Particles{CurrentChannel}(CurrentParticle).FrameApproved)),1)'*ErrorIntegral5,'.-','Color','blue');
%         p3 = plot(traceFigAxes,Frames(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
%             backGround3(Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.-','Color','blue');

try
    xlim(traceFigAxes,[min(traceFigTimeAxis),max(traceFigTimeAxis)]+[-1,1]);
catch
    %             error('Not sure what happened here. Problem with trace fig x lim. Talk to AR if you see this, please.');
end
% plotting anaphase boundaries ------------------------------------
% Section added by EL 10/11/18
currentYLimits = get(traceFigAxes,'YLim');  % Get the range of the y axis

% plotting all anaphase time points as vertical lines
for i = 1:length(anaphase)
    if ~lineFit
        currentAnaphaseBoundary = anaphase(i);
    else
        currentAnaphaseBoundary = anaphaseInMins(i) - priorAnaphaseInMins;
    end
    plot(traceFigAxes,ones(1,2).*currentAnaphaseBoundary,currentYLimits,...
        'LineWidth',2,'Color','black');
end
% End of anaphase boundary marking section
% -----------------------------------------------------------------

%----------------------------------------------------
% plotting prophase and metaphase boundaries
for i = 1:length(prophase)
    if ~lineFit
        currentProphase = prophase(i);
    else
        currentProphase = prophaseInMins(i) - priorAnaphaseInMins;
    end
    plot(traceFigAxes,ones(1,2).*currentProphase,currentYLimits,...
        'LineWidth',2,'Color','blue');
end

for i = 1:length(metaphase)
    if ~lineFit
        currentMetaphase = metaphase(i);
    else
        currentMetaphase = metaphaseInMins(i) - priorAnaphaseInMins;
    end
    plot(traceFigAxes,ones(1,2).*currentMetaphase,currentYLimits,...
        'LineWidth',2,'Color','yellow');
end

ylabel(traceFigAxes,'integrated intensity (a.u.)')
hold(traceFigAxes, 'off')
%         yyaxis(traceFigAxes,'right');
%         p3 = plot(traceFigAxes,Frames(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
%             backGround3(Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.-','Color','blue');
if plot3DGauss
    str1 = '3-slice';
    str2 = '3D-Gaussian fit';
else
    str1 = '1-slice';
    str2 = 'multi-slice';
end
if ~lineFit
    legend(traceFigAxes,[traceErrorBar1,traceErrorBar2],str1,str2)
    xlabel(traceFigAxes,'frame')
else
    legend(traceFigAxes,[traceErrorBar1,traceErrorBar2,fit1E],str1,str2,...
        ['fit slope: ', num2str(round(Coefficients(1))), ' a.u./min',newline,'time on: ',num2str(roots(Coefficients)), ' min'])
    xlabel(traceFigAxes,'time since anaphase (min)')
end


firstLine = ['Particle: ',num2str(CurrentParticle),'/',num2str(numParticles)];
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
title(traceFigAxes,axisTitle)
    
end

