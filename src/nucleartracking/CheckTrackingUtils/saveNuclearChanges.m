function saveNuclearChanges(cntState, DataFolder, FilePrefix, DropboxFolder)
    schnitzcells = cntState.schnitzcells;
    Ellipses = cntState.Ellipses;
    [schnitzcells, Ellipses] = correctSchnitzCellErrors(schnitzcells, Ellipses);
    FrameInfo = cntState.FrameInfo;

% If we only have one channel bring Particles back to the legacy format without any cells


save([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo', '-v6')



if cntState.UseHistoneOverlay
    if whos(var2str(schnitzcells)).bytes < 2E9
        save([DropboxFolder, filesep, FilePrefix(1:end-1), filesep, FilePrefix(1:end-1), '_lin.mat'], 'schnitzcells', '-v6')
    else
        save([DropboxFolder, filesep, FilePrefix(1:end-1), filesep, FilePrefix(1:end-1), '_lin.mat'], 'schnitzcells', '-v7.3', '-nocompression')
    end
end
save2([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,'Ellipses.mat'], Ellipses); 
disp('Schnitz Cells saved.')
end

