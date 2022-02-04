function CurrentFrame = previousFrameSegment(Particles, CurrentChannel, ...
    CurrentParticle, CurrentFrame)

AllFrames=Particles{CurrentChannel}(CurrentParticle).Frame;

%Frames in FrameRange that were not in this particle.
if CurrentFrame > min(AllFrames)
    
    if  CurrentFrame > min(AllFrames)
        RelevantFrames =[AllFrames(AllFrames < CurrentFrame) CurrentFrame];
        FrameSpaces=diff(RelevantFrames);
        if FrameSpaces(end)== 1
            RelFrameIndex=find(FrameSpaces>1,1, 'last');
            if ~isempty(RelFrameIndex)
                CurrentFrame=RelevantFrames(RelFrameIndex+1);
            else
                CurrentFrame=RelevantFrames(1);
            end
        else
            RelFrameIndex=find(FrameSpaces>1,2, 'last');
            if ~isempty(RelFrameIndex)
                CurrentFrame=RelevantFrames(RelFrameIndex(1)+1);
            else
                CurrentFrame=RelevantFrames(1);
            end
        end
    end
end

%end

%%

