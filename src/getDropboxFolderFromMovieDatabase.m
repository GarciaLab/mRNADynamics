function [dropboxFolderName, rowIndex] = getDropboxFolderFromMovieDatabase(movieDatabasePath, prefix, PREFIX_SEPARATOR, varargin)
  
  movieDatabase = csv2cell(movieDatabasePath, 'fromfile');
  
  if ~isempty(varargin)
    optionalResults = varargin{1};
    [dropboxFolderName, rowIndex] = getDropboxFolderFromMovieDatabaseContents(movieDatabase, prefix, PREFIX_SEPARATOR, optionalResults);
  else
    [dropboxFolderName, rowIndex] = getDropboxFolderFromMovieDatabaseContents(movieDatabase, prefix, PREFIX_SEPARATOR);
  end
  
end
