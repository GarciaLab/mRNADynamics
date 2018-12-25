function saveChanges(NChannels, Particles, Spots, SpotFilter, DataFolder, ...
    FrameInfo, UseHistoneOverlay, Threshold1, Threshold2, FilePrefix, ...
    schnitzcells, DropboxFolder)
%SAVECHANGES Summary of this function goes here
%   Detailed explanation goes here

%If we only have one channel bring Particles back to the legacy
%format without any cells
if NChannels==1
    Particles=Particles{1};
    Spots=Spots{1};
    SpotFilter=SpotFilter{1};
end

save([DataFolder,filesep,'FrameInfo.mat'],'FrameInfo')
if UseHistoneOverlay
    save([DataFolder,filesep,'Particles.mat'],'Particles','SpotFilter','Threshold1','Threshold2', '-v7.3')
    save([DataFolder,filesep,'Spots.mat'],'Spots', '-v7.3') %CS20170912 necessary for saving Spots.mat if >2GB
    save([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'],'schnitzcells', '-v7.3')
else
    save([DataFolder,filesep,'Particles.mat'],'Particles','SpotFilter','Threshold1','Threshold2', '-v7.3')
    save([DataFolder,filesep,'Spots.mat'],'Spots','-v7.3') %CS20170912 necessary for saving Spots.mat if >2GB
end
disp('Particles saved.')
if NChannels==1
    Particles={Particles};
    Spots = {Spots};
    SpotFilter = {SpotFilter};
end
end

