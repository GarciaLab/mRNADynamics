function [SourcePath, FISHPath, DropboxFolder, MS2CodePath, PreProcPath, configValues, movieDatabasePath] = DetermineLocalFolders(varargin)

  CONFIG_CSV_PATH = ['ComputerFolders.csv'];

  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');

  SourcePath = getConfigValue(configValues, 'SourcePath');
  FISHPath = getConfigValue(configValues, 'FISHPath');
  DropboxFolder = getConfigValue(configValues, 'DropboxFolder');
  MS2CodePath = getConfigValue(configValues, 'MS2CodePath');
  PreProcPath = getConfigValue(configValues, 'PreProcPath');
  movieDatabasePath = [DropboxFolder,'/MovieDatabase.csv'];

  if isempty(varargin) || isempty(varargin{1})
    warning('No Prefix specified. Using default Dropbox folder')
    return
  end
  
  PREFIX_SEPARATOR = '[\\\\/-]';
  PREFIX_REGEX = ['^.{10}', PREFIX_SEPARATOR, '.*$'];

  %% We need to look for the dropbox folder specified by the provided prefix
  Prefix = varargin{1};
  disp('Determining local folders using prefix:')
  disp(Prefix)
  if isempty(regexp(Prefix, PREFIX_REGEX))
    error('Prefix %s does not match "yyyy-mm-dd[/\\-]name". Please change it accordingly.', Prefix)
    % any 10 characters will work, not only yyyy-mm-dd,
    % but we enforce this in the error msg for simplicity
  end
  
  dropboxFolderName = getDropboxFolderFromMovieDatabase(movieDatabasePath, Prefix, PREFIX_SEPARATOR);

  %We can use the string "Default" or "DropboxFolder" for the default Dropbox folder
  if strcmpi(dropboxFolderName,'default')
    dropboxFolderName = 'DropboxFolder';
  end
  
  DropboxFolder = getConfigValue(configValues, dropboxFolderName);

end
