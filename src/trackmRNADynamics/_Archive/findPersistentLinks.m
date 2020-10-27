function [PrevIndices, NewIndices] = findPersistentLinks(ParticleStitchInfo,CurrentFrame)

PrevIndices = [];
NewIndices = [];
if ~isempty(ParticleStitchInfo)
  stitchFrameCell = ParticleStitchInfo.persistentLinkFrameCell;
  stitchIndexCell = ParticleStitchInfo.persistentLinkIndexCell;
  
  for s = 1:length(stitchFrameCell)
    joinFrames = stitchFrameCell{s};
    joinIndices = stitchIndexCell{s};
    if min(joinFrames)==CurrentFrame-1 && max(joinFrames)==CurrentFrame
      PrevIndices(end+1) = joinIndices(joinFrames==CurrentFrame-1);
      NewIndices(end+1) = joinIndices(joinFrames==CurrentFrame);
    end
  end
end