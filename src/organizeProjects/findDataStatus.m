function dataStatusFilename = findDataStatus(DropboxFolder)

dataStatusFilename = '';

dataStatusDir = dir([DropboxFolder,filesep,'Data*Status.*']);   %May or may not have space between Data and Status

if length(dataStatusDir) == 1
    dataStatusFilename = [DropboxFolder,filesep,dataStatusDir];
elseif length(dataStatusDir) > 1
    error(['More than one DataStatus.XLS found in folder ',DropboxFolder])
end