% For imaging sessions in which multiple embryos were imaged, separates out
% images into different folders.

% To do: create description of how to title/position files

function PrefixListSRCA = SeparateRawData

%Figure out the initial folders. We'll update the Drobpox one later on in the code.
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders; %#ok<NASGU,ASGLU>

OrigPathS = uigetdir(SourcePath, 'Select the folder (named after the imaged line) to be split into multiple folders');
SlashPositionR = find(OrigPathS == filesep);
DateFolderS = OrigPathS(SlashPositionR(end-1) + 1 : SlashPositionR(end) - 1);
LineFolderS = OrigPathS(SlashPositionR(end)+1 : end);

% Identifying unique prefixes of file names
DirD = dir([OrigPathS, filesep, '*.tif']);
FileListSRCA = {DirD.name};
NumFiles = length(FileListSRCA);
AllEmbryoPrefixesSRCA = cell(1, NumFiles);
for lFile = 1:NumFiles
    UnderscorePositionR = find(FileListSRCA{lFile} == '_');
    AllEmbryoPrefixesSRCA{lFile} = FileListSRCA{lFile}(1 : UnderscorePositionR(1)-1);
end
MoviePrefixSRCA = unique(AllEmbryoPrefixesSRCA);
NumMovies = length(MoviePrefixSRCA);

% Separate movie frames and low-zoom surface / midsagittal plane images
% into separate folders
cd(OrigPathS);
for lMovie = 1:NumMovies
    NewPathS = [OrigPathS, filesep, '..', filesep, LineFolderS, '_', MoviePrefixSRCA{lMovie}];
    if ~exist(NewPathS, 'dir')
        mkdir(NewPathS);
    end
    if ~exist([NewPathS, filesep, 'FullEmbryo'], 'dir')
        mkdir([NewPathS, filesep, 'FullEmbryo']);
    end
    SelectedFileD = dir([OrigPathS, filesep, MoviePrefixSRCA{lMovie}, '*.tif']);
    SelectedFileListSRCA = {SelectedFileD.name};
    NumSelectedFiles = length(SelectedFileListSRCA);
    for lSelectedFile = 1:NumSelectedFiles
        SelectedFileS = SelectedFileListSRCA{lSelectedFile};
        if strfind(SelectedFileS, 'Movie')
            movefile(SelectedFileS, NewPathS);
        else
            movefile(SelectedFileS, [NewPathS, filesep, 'FullEmbryo']);
        end
    end
end

% Delete original folder
if length(dir(OrigPathS)) == 2
    % This is a hack, because I don't think dir of a folder can be less
    % than 2 normally.
    cd('..');
    rmdir(OrigPathS)
else
    error('GregorLab:Error', 'Error: original folder has not been emptied.');
end

% Create prefix list
PrefixListSRCA = cell(1, NumMovies);
for lMovie = 1:NumMovies
    PrefixListSRCA{lMovie} = [DateFolderS, '-', LineFolderS, '_', MoviePrefixSRCA{lMovie}];
end

end