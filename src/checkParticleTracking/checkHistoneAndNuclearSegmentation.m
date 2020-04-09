%Check if we have the histone channel and we have done the nuclear segmentation.
function [Ellipses, UseHistoneOverlay, UseSchnitz] = checkHistoneAndNuclearSegmentation(PreProcPath, FilePrefix, NDigits, DropboxFolder, noHisOverlay, fish)

Ellipses = [];
UseHistoneOverlay = false;
Prefix = FilePrefix(1:end-1);

if exist([PreProcPath, filesep, FilePrefix(1:end - 1), filesep, ...
        FilePrefix(1:end - 1), '-His_', iIndex(1, NDigits), '.tif'], 'file') || exist([PreProcPath, filesep, Prefix, filesep, Prefix, '_hisMat.mat'], 'file')
    %(MT, 2018-02-11) Added support for lattice imaging with bad histone
    %channel, maybe temporary - FIX LATER
    if exist([DropboxFolder, filesep, FilePrefix(1:end - 1), filesep, 'Ellipses.mat'], 'file')
        load([DropboxFolder, filesep, FilePrefix(1:end - 1), filesep, 'Ellipses.mat'], 'Ellipses');
        UseHistoneOverlay = true;
    else
        warning('Ellipses.mat does not exist. Proceeding as though there is no Histone channel. If you expect a Histone channel, there is something wrong.')
        UseHistoneOverlay = false;
        Ellipses = [];
    end
end

%Check that we have the nuclear tracking done using schnitzcells
UseSchnitz = exist([DropboxFolder, filesep, FilePrefix(1:end - 1), filesep, FilePrefix(1:end - 1), '_lin.mat'], 'file');

if noHisOverlay || fish
    UseHistoneOverlay = false;
end

end