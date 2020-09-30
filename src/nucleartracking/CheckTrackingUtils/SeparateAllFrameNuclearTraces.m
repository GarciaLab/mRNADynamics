function schnitzcells=SeparateAllFrameNuclearTraces(CurrentNucleus, CurrentFrame, schnitzcells, FrameInfo, ncFrames)


%Separate the particle trace at the specified position
HoldCurrentNucleus = CurrentNucleus;
schnitzcells_inframe=[];
schnitzcells_inframe2=[];
numNuclei = length(schnitzcells);
for i=1:numNuclei
    if ~isempty(find(schnitzcells(i).frames ==CurrentFrame))
        schnitzcells_inframe= [schnitzcells_inframe,i];
    end
end
idx_adjustment = 0;
for i=1:length(schnitzcells_inframe)
    sc_idx = schnitzcells_inframe(i) + idx_adjustment;
    schnitzcells_inframe2(i) = sc_idx;
    frames_split1 = schnitzcells(sc_idx).frames(schnitzcells(sc_idx).frames < CurrentFrame);
    frames_split2 = schnitzcells(sc_idx).frames(schnitzcells(sc_idx).frames >= CurrentFrame);
    if (length(frames_split1) >= 1) & (length(frames_split2) >= 1)

       [CurrentNucleus,~, ~] = ...
            changeNucleus(sc_idx, schnitzcells, numNuclei);
       schnitzcells = SeparateNuclearTraces(CurrentNucleus, ...
        CurrentFrame, schnitzcells, FrameInfo, ncFrames); 
       idx_adjustment = idx_adjustment + 1;
    end
end


