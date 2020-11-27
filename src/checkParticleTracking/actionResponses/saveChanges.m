function saveChanges(NChannels, cptState, DataFolder, FilePrefix, DropboxFolder)
    Particles = cptState.Particles;
    Spots = cptState.Spots;
    SpotFilter = cptState.SpotFilter;
    FrameInfo = cptState.FrameInfo;
    schnitzcells = cptState.schnitzcells;

% If we only have one channel bring Particles back to the legacy format without any cells
if NChannels == 1
    Particles = Particles{1};
    Spots = Spots{1};
    SpotFilter = SpotFilter{1};
end

save([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo', '-v6')
save([DataFolder, filesep, 'Particles.mat'], 'Particles', 'SpotFilter', '-v6')

if whos(var2str(Spots)).bytes < 2E9
    save([DataFolder, filesep, 'Spots.mat'], 'Spots', '-v6')
else
    save([DataFolder, filesep, 'Spots.mat'], 'Spots', '-v7.3', '-nocompression')
end


if cptState.UseHistoneOverlay
    if whos(var2str(schnitzcells)).bytes < 2E9
        save([DropboxFolder, filesep, FilePrefix(1:end-1), filesep, FilePrefix(1:end-1), '_lin.mat'], 'schnitzcells', '-v6')
    else
        save([DropboxFolder, filesep, FilePrefix(1:end-1), filesep, FilePrefix(1:end-1), '_lin.mat'], 'schnitzcells', '-v7.3', '-nocompression')
    end
end

disp('Particles saved.')
end

