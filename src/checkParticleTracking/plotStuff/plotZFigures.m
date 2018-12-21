function [MaxZProfile, Frames] = plotZFigures(zProfileFigAxes, zTraceAxes, ExperimentType, ...
    xTrace, Spots, CurrentFrame, CurrentChannel, CurrentParticleIndex, ZSlices, ...
    CurrentZ, CurrentZIndex, PreviousParticle, CurrentParticle, ...
    PreviousChannel, Particles, MaxZProfile, Frames)
%PLOTZFIGURES Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(xTrace)

    %Get the z-DoG profile
    % check to see if Spots contains flag indicating type of
    % integration used
    title_string = '(max intensity)';
    %         if isfield(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex),'IntegralZ')
    if isfield(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex),'FixedAreaIntensity3')

        %             if Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).IntegralZ
        g = [-1 0 1];
        gaussFilter = exp(-g .^ 2 / (2 ));
        zprofinit = zeros(1, ZSlices);
        zprofinit(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).z) = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).FixedAreaIntensity;
        ZProfile= conv(gaussFilter,zprofinit);
        ZProfile = ZProfile(2:end-1);
        ZProfile = ZProfile(zprofinit~=0);
        title_string = '';
        IntegralZ_flag = 1;
        %             else
        %                 ZProfile=Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).CentralIntensity;
        %                 IntegralZ_flag = 0;
        %             end
    else
        ZProfile=Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).CentralIntensity;
        IntegralZ_flag = 0;
    end
    MaxZ=Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).brightestZ;

    plot(zProfileFigAxes, Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).z,...
        ZProfile,'.-k');

    hold(zProfileFigAxes,'on')
    if ~isempty(CurrentZIndex)
        plot(zProfileFigAxes,CurrentZ,ZProfile(CurrentZIndex),'ob')
    else
        plot(zProfileFigAxes,CurrentZ,CurrentZ,'or')
    end
    hold(zProfileFigAxes,'off')
    ylabel(zProfileFigAxes,'intensity(au)', 'FontSize',12);
    xlabel(zProfileFigAxes,'z-slice', 'FontSize',12);
    title(zProfileFigAxes,{'z-profile:';title_string},'FontSize',10)
end

%%%%%    BRIGHTEST Z-TRACE PLOT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmpi(ExperimentType,'inputoutput')
    %Only update the trace information if we have switched particles
    if (CurrentParticle~=PreviousParticle)||~exist('MaxZProfile', 'var')||(CurrentChannel~=PreviousChannel)|| CurrentFrame~=PreviousChannel
        PreviousParticle=CurrentParticle;
        Frames=PlotParticleTrace(CurrentParticle,Particles{CurrentChannel},Spots{CurrentChannel});
    end
    for  i = 1:length(Frames)
        MaxZProfile(i)=Spots{CurrentChannel}(Frames(i)).Fits...
            (Particles{CurrentChannel}(CurrentParticle).Index(i)).brightestZ;
    end
    plot(zTraceAxes,Frames(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
        MaxZProfile(Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.-k');
    hold(zTraceAxes, 'on')
    plot(zTraceAxes,Frames(Frames==CurrentFrame),MaxZProfile(Frames==CurrentFrame),'ob');
    hold(zTraceAxes, 'off')

    try
        xlim(zTraceAxes,[min(Frames)-1,max(Frames)+1]);
        ylim(zTraceAxes,[1,ZSlices+1])
    catch
        %             error('Not sure what happened here. Problem with trace fig x lim. Talk to AR if you see this, please.');
    end
    xlabel(zTraceAxes,'frame')
    ylabel(zTraceAxes,'z slice')
    title(zTraceAxes,'Brightest Z trace')
end
    
end

