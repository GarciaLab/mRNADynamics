function [Frames,Amp, PreviousParticle] = plotTraceInputOutput(traceFigAxes, ...
            FrameInfo, CurrentChannel, PreviousChannel, ...
    CurrentParticle, PreviousParticle, lastParticle, HideApprovedFlag, ...
    schnitzcells, Particles, numFrames, CurrentFrame, ZSlices, CurrentZ, Spots, ...
    Frames, Amp)
%PLOTTRACEINPUTOUTPUT Summary of this function goes here
%   Detailed explanation goes here

numParticles = length(Particles{CurrentChannel});

 %Only update the trace information if we have switched particles
if (CurrentParticle~=PreviousParticle)||~exist('Amp', 'var')||(CurrentChannel~=PreviousChannel) || lastParticle
    PreviousParticle=CurrentParticle;
    [Frames,Amp]=PlotParticleTrace(CurrentParticle,Particles{CurrentChannel},Spots{CurrentChannel});
end
cla(traceFigAxes, 'reset');
%we'll plot the spot intensity first on the left axis.
yyaxis(traceFigAxes,'left')
hold(traceFigAxes,'on')
plot(traceFigAxes,Frames(Particles{CurrentChannel}(CurrentParticle).FrameApproved),Amp(Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.-b')
xlabel(traceFigAxes,'frame')
ylabel(traceFigAxes,'transcript intensity (a.u.)')
plot(traceFigAxes,Frames(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),Amp(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.k')
plot(traceFigAxes,Frames(Frames==CurrentFrame),Amp(Frames==CurrentFrame),'ob')
% Should the nc boundary lines be added here?
% Search: plotting anaphase boundaries ------------------------------------
hold(traceFigAxes,'off')


%now we'll plot the input protein intensity on the right-hand axis.
yyaxis(traceFigAxes,'right')
plot(traceFigAxes,schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames,...
    max(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).Fluo,[],2),'r.-')
try
    xlim(traceFigAxes,[min(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames),max(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames)])
catch
    %             error('Not sure what happened here. Problem with trace fig x lim. Talk to AR if you see this, please.');
end
ylabel(traceFigAxes,'input protein intensity (a.u.)');
hold(traceFigAxes,'on')
plot(traceFigAxes,Frames(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),Amp(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.r')
hold(traceFigAxes,'off')

firstLine = ['Particle: ',num2str(CurrentParticle),'/',num2str(numParticles)];
secondLine = ['Frame: ',num2str(CurrentFrame),'/',num2str(numFrames),'    ',num2str(round(FrameInfo(CurrentFrame).Time)), 's'];
thirdLine = ['Z: ',num2str(CurrentZ),'/',num2str(ZSlices),', Ch: ',num2str(CurrentChannel)];

if isfield(FrameInfo, 'nc')
    FigureTitle={firstLine,...
        [secondLine,'    (nc',num2str(FrameInfo(CurrentFrame).nc),')'],...
        thirdLine};
else
    FigureTitle={firstLine,secondLine,thirdLine};
end

if HideApprovedFlag==1
    FigureTitle=[FigureTitle,', Showing non-flagged particles'];
elseif HideApprovedFlag==2
    FigureTitle=[FigureTitle,', Showing disapproved particles'];
end
title(traceFigAxes,FigureTitle)

end

