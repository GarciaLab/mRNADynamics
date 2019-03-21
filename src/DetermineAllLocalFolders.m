function [SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath, configValues, movieDatabasePath] = DetermineAllLocalFolders(varargin)

%Get default Dropbox folder location
CONFIG_CSV_PATH = 'ComputerFolders.csv';
configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');    
DefaultDropboxFolder = getConfigValue(configValues, 'DropboxFolder');
%Now get rest of folder locations
if(~isempty(varargin))   
    [SourcePath, FISHPath, DropboxFolder, MS2CodePath, PreProcPath, configValues, movieDatabasePath] = DetermineLocalFolders(varargin{1});
else
  	[SourcePath, FISHPath, DropboxFolder, MS2CodePath, PreProcPath, configValues, movieDatabasePath] = DetermineLocalFolders;
end
