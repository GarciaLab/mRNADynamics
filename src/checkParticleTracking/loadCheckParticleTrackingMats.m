% This function loads Particles.mat, Spots.mat and FrameInfo.mat into
% the corresponding workspace variables.
function [ParticleStitchInfo, Particles, SimParticles, SpotFilter, Spots,...
  FrameInfo, schnitzcells, Spots3D] = loadCheckParticleTrackingMats(DataFolder, PreProcPath, FilePrefix)

disp('Loading Particles.mat...')
load([DataFolder, filesep, 'Particles.mat'], 'Particles', 'SpotFilter')
load([DataFolder, filesep, 'ParticlesFull.mat'], 'SimParticles')
load([DataFolder, filesep, 'ParticleStitchInfo.mat'], 'ParticleStitchInfo')
disp('Particles.mat loaded')

schnitzcells = [];
schnitzPath = [DataFolder filesep, FilePrefix(1:end - 1), '_lin.mat'];
if  exist(schnitzPath, 'file')
    disp('Loading schnitzcells...')
    load(schnitzPath, 'schnitzcells');
    %Remove the schnitz fields that can give us problems potentially if
    %present. I don't know how this came to be, but it's for fields that
    %are not all that relevant. The fields are: approved, ang
    if isfield(schnitzcells, 'approved')
        schnitzcells = rmfield(schnitzcells, 'approved');
    end
    if isfield(schnitzcells, 'ang')
        schnitzcells = rmfield(schnitzcells, 'ang');
    end
    disp('schnitzcells loaded')
end


disp('Loading Spots.mat...')
if exist([DataFolder, filesep, 'SpotsMinimal.mat'], 'file')
    load([DataFolder, filesep, 'SpotsMinimal.mat'], 'SpotsMinimal');
    Spots = SpotsMinimal;
else
    load([DataFolder, filesep, 'Spots.mat'], 'Spots')
end
% NL: added for backwards compatibility
if ~iscell(Spots)
  Spots = {Spots};
end
if exist([DataFolder, filesep, 'Spots3D.mat'], 'file')
    load([DataFolder, filesep, 'Spots3D.mat'], 'Spots3D');
end

disp('Spots.mat loaded')
Spots3D = [];



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
