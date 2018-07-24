function [SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath, configValues, movieDatabasePath] = DetermineAllLocalFolders(varargin)

%Get the folders, including the default Dropbox one
[SourcePath, FISHPath, DefaultDropboxFolder, MS2CodePath, PreProcPath, configValues, movieDatabasePath] = DetermineLocalFolders;

%Now get the actual DropboxFolder according to the provided prefix
if(~isempty(varargin))
  [~, ~, DropboxFolder, ~, ~, ~, ~] = DetermineLocalFolders(varargin{1});
else
  DropboxFolder = DefaultDropboxFolder
end
