function dataStatusFolders = findDataStatus(dropboxFolders)

% dataStatusFolders = findDataStatus(DropboxFolders)
%
% DESCRIPTION
% Locates the DataStatus.xlsx file(s) inside the specified Dropbox (results)
% folder(s).
%
% INPUT
% dropboxFolders: Can be either a string containing the path to a single 
%                 DropboxFolder or a cell array containing the paths to 
%                 multiple DropboxFolders
%
% OPTIONS
% N/A
%
% OUTPUT
% dataStatusFolders: Cell array containing the DropboxFolder(s) where 
%                    DataStatus.xlsx file(s) were found
%
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 5/13/2020
% Origin: Functionalized from code originally included in LoadMS2Sets.m,
%         written by Hernan Garcia (hggarcia@berkeley.edu)
% Last Updated: N/A

%Initialize variables
dataStatusFilename = 'DataStatus.*';
dataStatusFolders = {};

%Allow for use with one or multiple DropboxFolders
if ischar(dropboxFolders)
    dropboxFolders = {dropboxFolders};
elseif ~iscell(dropboxFolders)
    error('DropboxFolders must be either a single string (char array) or a cell array of string(s) containing DropboxFolder path(s).')
end

%Search all DropboxFolders for DataStatus.xlsx files
for i = 1:length(dropboxFolders)
    currDropboxFolder = dropboxFolders{i};
    dataStatusDir = dir([currDropboxFolder,filesep,dataStatusFilename]);   
    
    %MT 5/15/20: Historically the code allowed either 'DataStatus.*' or 
    %'Data Status.*' (with a space between). I don't want to deal with 
    %that annoyance, so this is a check to enforce use of 'DataStatus.*' 
    %only
    dataSpaceStatusDir = dir([currDropboxFolder,filesep,'Data Status.*']);
    if ~isempty(dataSpaceStatusDir)
        error(['Found ''Data Status.xlsx'' (space between the words). Please rename to ''DataStatus.xlsx'' (no space) and re-try. Location: ', currDropboxFolder])
    end
    
    %Make sure we don't have multiple DataStatus.xlsx files in this 
    %DropboxFolder
    if length(dataStatusDir) == 1
        dataStatusFolders{end+1} = currDropboxFolder;
    elseif length(dataStatusDir) > 1
        error(['More than one DataStatus.xlsx found in folder ',currDropboxFolder, '. Please delete extra file(s) and re-try.'])
    end  
end

%Check that we found at least one DataStatus.xlsx file
if isempty(dataStatusFolders)
    fprintf(' %s\n,',dropboxFolders{:})
    error('No DataStatus.xlsx file in any of the locations listed above.')
end
