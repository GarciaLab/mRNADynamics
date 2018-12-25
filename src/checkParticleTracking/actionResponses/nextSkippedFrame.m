function CurrentFrame = nextSkippedFrame(Particles, CurrentChannel, ... 
    CurrentParticle, CurrentFrame)
%NEXTSKIPPEDFRAME Summary of this function goes here
%   Detailed explanation goes here

%This is the total frame range possible for this particle. Note
    %that we could still want to add spots at the beginning or end of
    %this range.
    FrameRange=Particles{CurrentChannel}(CurrentParticle).Frame(1):...
        Particles{CurrentChannel}(CurrentParticle).Frame(end);
    %Frames in FrameRange that were not in this particle.
    SkippedFrames=FrameRange(~ismember(FrameRange,...
        Particles{CurrentChannel}(CurrentParticle).Frame));
    %Find the next skipped frame and set it up
    CurrentFrame=min(SkippedFrames(SkippedFrames>CurrentFrame));

    %If there is no next empty frame, jump to the last frame of this
    %particle
    if isempty(CurrentFrame)
        CurrentFrame=Particles{CurrentChannel}(CurrentParticle).Frame(end);
    end
end

