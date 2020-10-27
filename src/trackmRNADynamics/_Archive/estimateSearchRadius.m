function searchRadius = estimateSearchRadius(Spots,ncVec,maxSearchRadius,PixelSize,matchFraction,Channel)
  % get number of frames
  nFrames = length(Spots{Channel});
  % estimate appropriate jump threshold
  nSpotsVec = NaN(1,length(Spots{Channel}));
  rng(234);
  for s = 1:length(Spots{Channel})
    nSpotsVec(s) = length(Spots{Channel}(s).Fits);
  end
  % select frames to use as test cases
  testFrames = randsample(1:length(Spots{Channel}),ceil(nFrames/33),true,nSpotsVec.*[diff(ncVec)==0 0]);
  % inialize vector to track nearest neighbor distances
  nnVec = [];
  % itereate through frames
  for testFrame = testFrames
    % Get spot positions
    [PrevSpotsX, PrevSpotsY, ~] = SpotsXYZ(Spots{Channel}(testFrame));
    [NewSpotsX, NewSpotsY, ~] = SpotsXYZ(Spots{Channel}(testFrame+1));
    nnVec = [nnVec min(sqrt((NewSpotsX'-PrevSpotsX).^2 + (NewSpotsY'- PrevSpotsY).^2)*PixelSize)];
  end
  % calculate search radius 
  searchRadius = min([maxSearchRadius, prctile(nnVec,matchFraction)]);