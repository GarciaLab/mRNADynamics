% This function loads Particles.mat, Spots.mat and FrameInfo.mat into
% the corresponding workspace variables
function [Particles, SpotFilter, Spots, FrameInfo] = loadCheckParticleTrackingMats(DataFolder, PreProcPath)

disp('Loading Particles.mat...')
  load([DataFolder, filesep, 'Particles.mat'], 'Particles', 'SpotFilter')
  disp('Particles.mat loaded')
  disp('Loading Spots.mat...')
  load([DataFolder, filesep, 'Spots.mat'], 'Spots')
  disp('Spots.mat loaded')

  %Check that FrameInfo exists
  if exist([DataFolder, filesep, 'FrameInfo.mat'], 'file')
    load([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo')
  else
    warning('No FrameInfo.mat found. Trying to continue')
    %Adding frame information
    DHis = dir([PreProcPath, filesep, FilePrefix(1:end - 1), filesep, '*His*.tif']);
    FrameInfo(length(DHis)).nc = [];
    %Adding information
    Dz = dir([PreProcPath, filesep, FilePrefix(1:end - 1), filesep, FilePrefix(1:end - 1), '*001*.tif']);
    NumberSlices = length(Dz) - 1;

    for i = 1:length(FrameInfo)
      FrameInfo(i).NumberSlices = NumberSlices;
    end

  end
end