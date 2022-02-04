function CurrentFrame = nextFrameSegment(Particles, CurrentChannel, ...
    CurrentParticle, CurrentFrame)

AllFrames=Particles{CurrentChannel}(CurrentParticle).Frame;
%Frames in FrameRange that were not in this particle.
if CurrentFrame < max(AllFrames)
    
    RelevantFrames =[CurrentFrame AllFrames(AllFrames > CurrentFrame)];
    
    %Find the next chunk of continuous frame and set it up
    FrameSpaces=diff(RelevantFrames);
    RelFrameIndex=find(FrameSpaces>1,1)+1;
    if ~isempty(RelFrameIndex)
        CurrentFrame=RelevantFrames(RelFrameIndex);
    end
end

end


