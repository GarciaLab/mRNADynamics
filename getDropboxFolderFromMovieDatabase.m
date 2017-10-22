function [dropboxFolderName, rowIndex] = getDropboxFolderFromMovieDatabase(movieDatabasePath, prefix, PREFIX_SEPARATOR)
  movieDatabase = csv2cell(movieDatabasePath, 'fromfile');
  [dropboxFolderName, rowIndex] = getDropboxFolderFromMovieDatabaseContents(movieDatabase, prefix, PREFIX_SEPARATOR);
end
