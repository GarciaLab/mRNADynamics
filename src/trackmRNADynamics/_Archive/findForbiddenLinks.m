function [PrevIndices, NewIndices] = findForbiddenLinks(ParticleStitchInfo,CurrentFrame)

PrevIndices = [];
NewIndices = [];
if ~isempty(ParticleStitchInfo)
  stitchFrameCell = ParticleStitchInfo.forbiddenLinkFrameCell;
  stitchIndexCell = ParticleStitchInfo.forbiddenLinkIndexCell;
  
  for s = 1:length(stitchFrameCell)
    splitFrames = stitchFrameCell{s};
    splitIndices = stitchIndexCell{s};
    if min(splitFrames)==CurrentFrame-1 && max(splitFrames)==CurrentFrame
      PrevIndices(end+1) = splitIndices(splitFrames==CurrentFrame-1);
      NewIndices(end+1) = splitIndices(splitFrames==CurrentFrame);
    end
  end
end