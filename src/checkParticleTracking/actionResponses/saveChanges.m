function saveChanges(NChannels, Particles, Spots, SpotFilter, DataFolder, ...
    FrameInfo, UseHistoneOverlay, FilePrefix, ...
    schnitzcells, DropboxFolder)

    % If we only have one channel bring Particles back to the legacy format without any cells
    if NChannels == 1
        Particles = Particles{1};
        Spots = Spots{1};
        SpotFilter = SpotFilter{1};
    end

    save([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo')
    
    % CS20170912, saves as 7.3 version, necessary for saving mats if >2GB
    if UseHistoneOverlay
        save([DataFolder, filesep, 'Particles.mat'], 'Particles', 'SpotFilter', '-v7.3')
        save([DataFolder, filesep, 'Spots.mat'], 'Spots', '-v7.3') 
        save([DropboxFolder, filesep, FilePrefix(1:end-1), filesep, FilePrefix(1:end-1), '_lin.mat'], 'schnitzcells', '-v7.3')
    else
        save([DataFolder, filesep, 'Particles.mat'], 'Particles', 'SpotFilter', '-v7.3')
        save([DataFolder, filesep, 'Spots.mat'], 'Spots', '-v7.3')
    end
    
    disp('Particles saved.')
end

