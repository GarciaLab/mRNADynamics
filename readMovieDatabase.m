function [SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
    Folder, Prefix, ExperimentType,Channel1,Channel2,OutputFolder]...
= readMovieDatabase(PrefixOverrideFlag)

%Figure out the initial folders. We'll update the Drobpox one later on in the code.
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath]=...
    DetermineLocalFolders;

%Get the folder with the data
if ~PrefixOverrideFlag
    Folder=uigetdir(SourcePath,'Select folder with data');
else 
    Folder = '';
end

%Get the information from the last two folders in the structure
if ~PrefixOverrideFlag
    SlashPositions=strfind(Folder,filesep);
    
    Prefix=[Folder((SlashPositions(end-1)+1):(SlashPositions(end)-1)),'-',...
        Folder((SlashPositions(end)+1):(end))];
else 
    Prefix = PrefixOverrideFlag;
end

%What type of experiment are we dealing with? Get this out of
%MovieDatabase.xlsx
[XLSNum,XLSTxt]=xlsread([DropboxFolder,filesep,'MovieDatabase.xlsx']);
ExperimentTypeColumn=find(strcmp(XLSTxt(1,:),'ExperimentType'));
DataFolderColumn=find(strcmp(XLSTxt(1,:),'DataFolder'));
Channel1Column=find(strcmp(XLSTxt(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSTxt(1,:),'Channel2'));


Dashes=findstr(Prefix,'-');
PrefixRow=find(strcmp(XLSTxt(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
if isempty(PrefixRow)
    PrefixRow=find(strcmp(XLSTxt(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)
        error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
    end
end

ExperimentType=XLSTxt(PrefixRow,ExperimentTypeColumn);
Channel1=XLSTxt(PrefixRow,Channel1Column);
Channel2=XLSTxt(PrefixRow,Channel2Column);

%Set the destination folders
OutputFolder=[DropboxFolder,filesep,Prefix];

[~,~,DropboxFolder,~,~]=...
    DetermineLocalFolders(Prefix);