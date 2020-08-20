function [rootFolderName, rowIndex] = getRootFolderFromMovieDatabaseContents(movieDatabase, prefix, PREFIX_SEPARATOR)
    movieDatabaseHeaderRow = movieDatabase(1, :);

    dataFolderColumnIndex = findColumnIndex(movieDatabaseHeaderRow, 'DataFolder');
    dataFolderColumn = movieDatabase(:, dataFolderColumnIndex);

    rootFolderColumnIndex = findColumnIndex(movieDatabaseHeaderRow, 'RootFolder');
    if isempty(rootFolderColumnIndex)
        rootFolderName = 'noFolder';
        rowIndex = NaN;
    else
        namestart = find(isletter(prefix));
        if (isempty(namestart))
          % If the prefix name does not contain any letters, defaults to 12 
          % which is the usual size after the date.
          namestart = 12;
        end

        namestart = namestart(1); %Index of first letter in prefix name, i.e. start of dataset name

        dash_indices = find(prefix(1:namestart)=='-'); %Get indices of dashes before the start of prefix name
        sep_dash = dash_indices(end); %Get index of dash serving as prefix separator

        indexArray = regexpi(dataFolderColumn, ['^', prefix(1:(sep_dash-1)), PREFIX_SEPARATOR, strrep(prefix((sep_dash+1):end), '+', '\+'), '$']);
        rowIndex = find(not(cellfun('isempty', indexArray)));

        rootFolderNameCell = movieDatabase(rowIndex, rootFolderColumnIndex);

        if isempty(rootFolderNameCell)
        error(['Data set "', prefix, '" not found in MovieDatabase.csv'])
        end

        rootFolderName = rootFolderNameCell{1};
    end
end
