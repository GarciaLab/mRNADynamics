function Particles = detrendZ(Particles,FrameInfo)

if isfield(FrameInfo,'zPosition')      
  % need to make sure this field is now a permanent feature
  zPosVec = [FrameInfo.zPosition]*1e6 / FrameInfo(1).ZStep;
  zPosVec = zPosVec - zPosVec(1);
  frameIndex = 1:length(zPosVec);
  % generate new det-trended z variable
  for p = 1:length(Particles)
    fVec = Particles(p).Frame;
    Particles(p).zPosDetrended = Particles(p).zPos - zPosVec(ismember(frameIndex,fVec));
  end
else
  % get list of all frames and corresponding z positions
  frameVec = [Particles.Frame];
  zPosVec = [Particles.zPos];

  % get iteratable frame  list
  frameIndex = unique(frameVec);
  avgZProfile = NaN(size(frameIndex));

  for i = 1:length(frameIndex)
    avgZProfile(i) = nanmean(zPosVec(frameVec==frameIndex(i)));
  end

  % generate new det-trended z variable
  for p = 1:length(Particles)
    fVec = Particles(p).Frame;
    Particles(p).zPosDetrended = Particles(p).zPos - avgZProfile(ismember(frameIndex,fVec));
  end

end