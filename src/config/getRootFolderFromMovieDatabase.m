function [RootFolderName, rowIndex] = getRootFolderFromMovieDatabase(movieDatabase, prefix, PREFIX_SEPARATOR)
  
if ischar(movieDatabase) %accept the input as either the database itself or a path to the database
    movieDatabase = csv2cell(movieDatabasePath, 'fromfile');
  end
  
  [RootFolderName, rowIndex] = getRootFolderFromMovieDatabaseContents(movieDatabase, prefix, PREFIX_SEPARATOR);
end