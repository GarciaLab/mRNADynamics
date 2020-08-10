function [PrevIndices, NewIndices] = findPersistentLinks(ParticleStitchInfo,CurrentFrame)

PrevIndices = [];
NewIndices = [];
if ~isempty(ParticleStitchInfo)
  stitchFrameCell = {{61,62} {55,57}, {61, 65}};%[ParticleStitchInfo{1}.persistentLinkFrames];
  stitchIndexCell = [ParticleStitchInfo{1}.persistentLinkIndices];
  
  for s = 1:length(stitchFrameCell)
    joinFrames = [stitchFrameCell{s}{:}];
    joinIndices = [stitchindexCell{s}{:}];
    if min(joinFrames)==CurrentFrame-1 && max(joinFrames)==CurrentFrame
      PrevIndices(end+1) = joinIndices(joinFrames==CurrentFrame-1);
      NewIndices(end+1) = joinIndices(joinFrames==CurrentFrame);
    end
  end
end