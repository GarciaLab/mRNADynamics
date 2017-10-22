function [dropboxFolderName, rowIndex] = getDropboxFolderFromMovieDatabaseContents(movieDatabase, prefix, PREFIX_SEPARATOR)
  movieDatabaseHeaderRow = movieDatabase(1, :);

  dataFolderColumnIndex = findColumnIndex(movieDatabaseHeaderRow, 'DataFolder');
  dataFolderColumn = movieDatabase(:, dataFolderColumnIndex);

  dropboxFolderColumnIndex = findColumnIndex(movieDatabaseHeaderRow, 'DropboxFolder');

  indexArray = regexpi(dataFolderColumn, ['^', prefix(1:10), PREFIX_SEPARATOR, prefix(12:end), '$']);
  rowIndex = find(not(cellfun('isempty', indexArray)));

  dropboxFolderNameCell = movieDatabase(rowIndex, dropboxFolderColumnIndex);

  if isempty(dropboxFolderNameCell)
    error(['Data set "', prefix, '" not found in MovieDatabase.csv'])
  end

  dropboxFolderName = dropboxFolderNameCell{1};
end
