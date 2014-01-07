% For imaging sessions in which multiple embryos were imaged, separates out
% images into different folders.

% To do: create description of how to title/position files

function PrefixListRCA = SeparateRawData

%Figure out the initial folders. We'll update the Drobpox one later on in the code.
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders;

OrigPathS = uigetdir(SourcePath, 'Select the folder (named after the imaged line) to be split into multiple folders');
SlashPositionR = find(OrigPathS == filesep);
LineFolderS = OrigPathS(SlashPositionR(end)+1 : end);
%DateFolderS = OrigPathS(SlashPositionR(end-1) + 1 : SlashPositionR(end) - 1);

% Identifying unique prefixes of file names
DirD = dir([OrigPathS, filesep, '*.tif']);
FileListSRCA = {DirD.name};
NumFiles = length(FileListSRCA);
AllEmbryoPrefixesSRCA = cell(1, NumFiles);
for lFile = 1:NumFiles
    UnderscorePositionR = find(FileListSRCA{lFile} == '_');
    AllEmbryoPrefixesSRCA{lFile} = FileListSRCA{lFile}(1 : UnderscorePositionR(1)-1);
end
MoviePrefixSRCA = unique(AllEmbryoPrefixesSRCA{lFile});
NumMovies = length(MoviePrefixSRCA);

% Separate movie frames and low-zoom surface / midsagittal plane images
% into separate folders
for lMovie = 1:NumMovies
    
end

end