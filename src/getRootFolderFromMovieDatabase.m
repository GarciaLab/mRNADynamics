function [RootFolderName, rowIndex] = getRootFolderFromMovieDatabase(movieDatabasePath, prefix, PREFIX_SEPARATOR)
  movieDatabase = csv2cell(movieDatabasePath, 'fromfile');
  [RootFolderName, rowIndex] = getRootFolderFromMovieDatabaseContents(movieDatabase, prefix, PREFIX_SEPARATOR);
end