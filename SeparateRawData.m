% For imaging sessions in which multiple embryos were imaged, separates out
% images into different folders.

% To do: create description of how to title/position files

function PrefixListRCA = SeparateRawData

%Figure out the initial folders. We'll update the Drobpox one later on in the code.
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders;

OrigPathS = uigetdir(SourcePath, 'Select folder with data to be split into multiple folders');
DirD = dir([OrigPathS, filesep, '*.tif']);
FileListRCA = {DirD.name};
NumFiles = length(FileListRCA);

for lFile = 1:NumFiles
    
end

end