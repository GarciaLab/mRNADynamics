function [SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
    Folder, Prefix, ExperimentType,Channel1,Channel2,OutputFolder]...
= readMovieDatabase(PrefixOverrideFlag)

%Figure out the initial folders. We'll update the Drobpox one later on in the code.
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath, configValues, movieDatabasePath]=...
    DetermineLocalFolders;

%Get the folder with the data
if ~PrefixOverrideFlag
    Folder = uigetdir(SourcePath,'Select folder with data');

    %Get the information from the last two folders in the structure
    SlashPositions = strfind(Folder,filesep);
    Prefix = [Folder((SlashPositions(end-1)+1):(SlashPositions(end)-1)),'-',...
        Folder((SlashPositions(end)+1):(end))];
else 
    Folder = '';
    Prefix = PrefixOverrideFlag;
end


%What type of experiment are we dealing with? Get this out of MovieDatabase
movieDatabase = csv2cell(movieDatabasePath, 'fromfile');
movieDatabaseHeaderRow = movieDatabase(1, :);
ExperimentTypeColumn = findColumnIndex(movieDatabaseHeaderRow, 'ExperimentType')
Channel1Column = findColumnIndex(movieDatabaseHeaderRow, 'Channel1')
Channel2Column = findColumnIndex(movieDatabaseHeaderRow, 'Channel2')

[DropboxFolder, PrefixRow] = getDropboxFolderFromMovieDatabase(movieDatabasePath, Prefix, '[\\\\/-]')

ExperimentType = movieDatabase(PrefixRow, ExperimentTypeColumn);
%ExperimentType = ExperimentType{1}
Channel1 = movieDatabase(PrefixRow, Channel1Column);
%Channel1 = Channel1{1}
Channel2 = movieDatabase(PrefixRow, Channel2Column);
%Channel2 = Channel2{1}

[~,~,DropboxFolder,~,~] = DetermineLocalFolders(Prefix);

%Set the destination folders
OutputFolder = [DropboxFolder, filesep, Prefix]
