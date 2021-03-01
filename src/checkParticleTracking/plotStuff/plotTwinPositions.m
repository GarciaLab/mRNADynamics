function plotTwinPositions(zProfileFigAxes, zTraceAxes, cptState)

MinTwinX = NaN;
MaxTwinX = NaN;
MinTwinY = NaN;
MaxTwinY = NaN;
MinTwinZ = NaN;
MaxTwinZ = NaN;
MinTwinFrame = NaN;
MaxTwinFrame = NaN;


twinXhandle = zProfileFigAxes{1}.Children(end-1);
twinYhandle = zProfileFigAxes{2}.Children(end-1);
twinZhandle = zTraceAxes.Children(end-1);
twinXPointHandle = zProfileFigAxes{1}.Children(end-3);
twinYPointHandle = zProfileFigAxes{2}.Children(end-3);
twinZPointHandle = zTraceAxes.Children(end-3);

if cptState.UseTwinTraces & ~isempty(cptState.TwinParticle)
    TwinXprofile = cptState.Particles{cptState.CurrentChannelIndex}(cptState.TwinParticle).xPos;
    TwinYprofile = cptState.Particles{cptState.CurrentChannelIndex}(cptState.TwinParticle).yPos;
    TwinZprofile = cptState.Particles{cptState.CurrentChannelIndex}(cptState.TwinParticle).zPos;
    ApprovedTwinFrames = cptState.Particles{cptState.CurrentChannelIndex}(cptState.TwinParticle).FrameApproved;
    TwinFrames = cptState.Particles{cptState.CurrentChannelIndex}(cptState.TwinParticle).Frame;
    if ~isempty(TwinXprofile)
        MinTwinX = min(TwinXprofile);
        MaxTwinX = max(TwinXprofile);
        MinTwinY = min(TwinYprofile);
        MaxTwinY = max(TwinYprofile);
        MinTwinZ = min(TwinZprofile);
        MaxTwinZ = max(TwinZprofile);
        MinTwinFrame = min(TwinFrames);
        MaxTwinFrame = max(TwinFrames);
    end
    
    

    
    set(twinXhandle, 'XData', TwinFrames(ApprovedTwinFrames),...
        'YData', TwinXprofile(ApprovedTwinFrames), 'HandleVisibility', 'on', 'Visible', 'on');
    set(twinYhandle, 'XData', TwinFrames(ApprovedTwinFrames),...
        'YData', TwinYprofile(ApprovedTwinFrames), 'HandleVisibility', 'on', 'Visible', 'on');
    set(twinZhandle, 'XData', TwinFrames(ApprovedTwinFrames),...
        'YData', TwinZprofile(ApprovedTwinFrames), 'HandleVisibility', 'on', 'Visible', 'on');
    if ~isempty(find(TwinFrames == cptState.CurrentFrame, 1))
        set(twinXPointHandle, 'XData', TwinFrames(TwinFrames == cptState.CurrentFrame),...
            'YData', TwinXprofile(TwinFrames == cptState.CurrentFrame), 'HandleVisibility', 'on', 'Visible', 'on');
        set(twinYPointHandle, 'XData', TwinFrames(TwinFrames == cptState.CurrentFrame),...
            'YData', TwinYprofile(TwinFrames == cptState.CurrentFrame), 'HandleVisibility', 'on', 'Visible', 'on');
        set(twinZPointHandle, 'XData', TwinFrames(TwinFrames == cptState.CurrentFrame),...
            'YData', TwinZprofile(TwinFrames == cptState.CurrentFrame), 'HandleVisibility', 'on', 'Visible', 'on');
    else
        set(twinXPointHandle, 'Visible', 'off');
        set(twinYPointHandle, 'Visible', 'off');
        set(twinZPointHandle, 'Visible', 'off');
    end
elseif cptState.UseTwinTraces  & isempty(cptState.TwinParticle)
    set(twinXPointHandle, 'Visible', 'off');
    set(twinYPointHandle, 'Visible', 'off');
    set(twinZPointHandle, 'Visible', 'off');
    set(twinXhandle, 'Visible', 'off');
    set(twinYhandle, 'Visible', 'off');
    set(twinZhandle, 'Visible', 'off');
end

CurrentXprofile = cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).xPos;
CurrentYprofile = cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).yPos;
CurrentZprofile = cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).zPos;
ApprovedFrames = cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).FrameApproved;
CurrentFrames = cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).Frame;

MinX = min(CurrentXprofile);
MaxX = max(CurrentXprofile);
MinY = min(CurrentYprofile);
MaxY = max(CurrentYprofile);
MinZ = min(CurrentZprofile);
MaxZ = max(CurrentZprofile);
MinFrame = min(CurrentFrames);
MaxFrame = max(CurrentFrames);



currentXhandle = zProfileFigAxes{1}.Children(end);
currentYhandle = zProfileFigAxes{2}.Children(end);
currentZhandle = zTraceAxes.Children(end);
currentXPointHandle = zProfileFigAxes{1}.Children(end-2);
currentYPointHandle = zProfileFigAxes{2}.Children(end-2);
currentZPointHandle = zTraceAxes.Children(end-2);

set(currentXhandle, 'XData', CurrentFrames(ApprovedFrames),...
    'YData', CurrentXprofile(ApprovedFrames), 'HandleVisibility', 'on', 'Visible', 'on');
set(currentYhandle, 'XData', CurrentFrames(ApprovedFrames),...
    'YData', CurrentYprofile(ApprovedFrames), 'HandleVisibility', 'on', 'Visible', 'on');
set(currentZhandle, 'XData', CurrentFrames(ApprovedFrames),...
    'YData', CurrentZprofile(ApprovedFrames), 'HandleVisibility', 'on', 'Visible', 'on');
if ~isempty(find(CurrentFrames == cptState.CurrentFrame, 1))
    set(currentXPointHandle, 'XData', CurrentFrames(CurrentFrames == cptState.CurrentFrame),...
        'YData', CurrentXprofile(CurrentFrames == cptState.CurrentFrame), 'HandleVisibility', 'on', 'Visible', 'on');
    set(currentYPointHandle, 'XData', CurrentFrames(CurrentFrames == cptState.CurrentFrame),...
        'YData', CurrentYprofile(CurrentFrames == cptState.CurrentFrame), 'HandleVisibility', 'on', 'Visible', 'on');
    set(currentZPointHandle, 'XData', CurrentFrames(CurrentFrames == cptState.CurrentFrame),...
        'YData', CurrentZprofile(CurrentFrames == cptState.CurrentFrame), 'HandleVisibility', 'on', 'Visible', 'on');
else
    set(currentXPointHandle, 'Visible', 'off');
    set(currentYPointHandle, 'Visible', 'off');
    set(currentZPointHandle, 'Visible', 'off');
    set(twinXPointHandle, 'Visible', 'off');
    set(twinYPointHandle, 'Visible', 'off');
    set(twinZPointHandle, 'Visible', 'off');
end


