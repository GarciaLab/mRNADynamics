%Check if we have the histone channel and we have done the nuclear segmentation.
function [Ellipses, UseHistoneOverlay, UseSchnitz] = checkHistoneAndNuclearSegmentation(...
    PreProcPath, FilePrefix, NDigits, DropboxFolder, noHisOverlay, fish)

Ellipses = [];
UseHistoneOverlay = false;
Prefix = FilePrefix(1:end-1);

liveExperiment = LiveExperiment(Prefix);

try
    Ellipses = getEllipses(liveExperiment);
    UseHistoneOverlay = true;
catch
    Ellipses = [];
    UseHistoneOverlay = false;
end

%Check that we have the nuclear tracking done using schnitzcells
UseSchnitz = exist([DropboxFolder, filesep,...
    FilePrefix(1:end - 1), filesep, FilePrefix(1:end - 1), '_lin.mat'], 'file');

if noHisOverlay || fish
    UseHistoneOverlay = false;
end

end