function CurrentFrame = previousNuclearUnapprovedFrame(schnitzcells,  ... 
    CurrentNucleus, CurrentFrame)
%PREVIOUSSKIPPEDFRAME Summary of this function goes here
%   Detailed explanation goes here
%This is the total frame range possible for this particle. Note
%that we could still want to add spots at the beginning or end of
%this range.
OldFrame = CurrntFrame;
FrameRange=schnitzcells(CurrentNucleus).frames(1):...
    schnitzcells(CurrentNucleus).frames(end);
%Frames in FrameRange that were not in this particle.
SkippedFrames=FrameRange(~ismember(FrameRange,...
    schnitzcells(CurrentNucleus).frames(schnitzcells(CurrentNucleus).FrameApproved)));
%Find the next skipped frame and set it up
CurrentFrame=max(SkippedFrames(SkippedFrames<CurrentFrame));

%If there is no next empty frame, jump to the last frame of this
%particle
if isempty(CurrentFrame)
    CurrentFrame=OldFrame;
end
end

