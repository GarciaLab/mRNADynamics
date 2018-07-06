function [dropboxFolderName, rowIndex] = getDropboxFolderFromMovieDatabaseContents(movieDatabase, prefix, PREFIX_SEPARATOR)
  movieDatabaseHeaderRow = movieDatabase(1, :);

  dataFolderColumnIndex = findColumnIndex(movieDatabaseHeaderRow, 'DataFolder');
  dataFolderColumn = movieDatabase(:, dataFolderColumnIndex);

  dropboxFolderColumnIndex = findColumnIndex(movieDatabaseHeaderRow, 'DropboxFolder');
  
%   namestart = find(isletter(prefix));
%   namestart = namestart(1); %Index of first letter in prefix name, i.e. start of dataset name
% 
%   indexArray = regexpi(dataFolderColumn, ['^', prefix(1:(namestart-2)), PREFIX_SEPARATOR, strrep(prefix(namestart:end), '+', '\+'), '$']);
    indexArray = regexpi(dataFolderColumn, ['^', prefix(1:10), PREFIX_SEPARATOR, strrep(prefix(12:end),'+', '\+'), '$']);  
    rowIndex = find(not(cellfun('isempty', indexArray)));

  dropboxFolderNameCell = movieDatabase(rowIndex, dropboxFolderColumnIndex);

  if isempty(dropboxFolderNameCell)
    error(['Data set "', prefix, '" not found in MovieDatabase.csv'])
  end

  dropboxFolderName = dropboxFolderNameCell{1};
end
