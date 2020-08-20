function [dropboxFolderName, rowIndex, movieDatabase] = getDropboxFolderFromMovieDatabase(movieDatabase, prefix, PREFIX_SEPARATOR, varargin)
  
if ischar(movieDatabase) %accept the input as either the database itself or a path to the database
  movieDatabase = csv2cell(movieDatabase, 'fromfile');
end

  if ~isempty(varargin)
    optionalResults = varargin{1};
  else
      optionalResults = '';
  end
  
  [dropboxFolderName, rowIndex] = getDropboxFolderFromMovieDatabaseContents(movieDatabase, prefix, PREFIX_SEPARATOR, optionalResults);

  
end
