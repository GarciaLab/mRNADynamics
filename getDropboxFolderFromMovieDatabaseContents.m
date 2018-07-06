function [dropboxFolderName, rowIndex] = getDropboxFolderFromMovieDatabaseContents(movieDatabase, prefix, PREFIX_SEPARATOR)
  movieDatabaseHeaderRow = movieDatabase(1, :);

  dataFolderColumnIndex = findColumnIndex(movieDatabaseHeaderRow, 'DataFolder');
  dataFolderColumn = movieDatabase(:, dataFolderColumnIndex);

  dropboxFolderColumnIndex = findColumnIndex(movieDatabaseHeaderRow, 'DropboxFolder');
  
  namestart = find(isletter(prefix));
  namestart = namestart(1); %Index of first letter in prefix name, i.e. start of dataset name
  
  dash_indices = find(prefix(1:namestart)=='-'); %Get indices of dashes before the start of prefix name
  sep_dash = dash_indices(end); %Get index of dash serving as prefix separator
 
  indexArray = regexpi(dataFolderColumn, ['^', prefix(1:(sep_dash-1)), PREFIX_SEPARATOR, strrep(prefix((sep_dash+1):end), '+', '\+'), '$']);
  rowIndex = find(not(cellfun('isempty', indexArray)));

  dropboxFolderNameCell = movieDatabase(rowIndex, dropboxFolderColumnIndex);

  if isempty(dropboxFolderNameCell)
    error(['Data set "', prefix, '" not found in MovieDatabase.csv'])
  end

  dropboxFolderName = dropboxFolderNameCell{1};
end
