function [dropboxFolderName, rowIndex] = getDropboxFolderFromMovieDatabaseContents(movieDatabase, prefix, PREFIX_SEPARATOR, varargin)
  
  if ~isempty(varargin)
     optionalResults = varargin{1};
  else
      optionalResults = '';
  end

  movieDatabaseHeaderRow = movieDatabase(1, :);

  dataFolderColumnIndex = findColumnIndex(movieDatabaseHeaderRow, 'DataFolder');
  dataFolderColumn = movieDatabase(:, dataFolderColumnIndex);

  dropboxFolderColumnIndex = findColumnIndex(movieDatabaseHeaderRow, 'DropboxFolder');
  
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

  dropboxFolderNameCell = movieDatabase(rowIndex, dropboxFolderColumnIndex);

  if isempty(dropboxFolderNameCell)
    error(['Data set "', prefix, '" not found in MovieDatabase.csv'])
  end

  if length(dropboxFolderNameCell) > 1
      if ~isempty(optionalResults)
          dropboxFolderName = optionalResults;
          rowIndex = rowIndex(2);
      else
          dropboxFolderName = dropboxFolderNameCell{1};
          rowIndex = rowIndex(1);
      end 
  else
     dropboxFolderName = dropboxFolderNameCell{1};
     rowIndex = rowIndex(1);
  end
     
%      f = figure();
%      c = uicontrol(f,'Style','popupmenu');
%      c.Position = [20 75 60 20];
%      c.String = {'dummy'};
%      c.Callback = @selection;
%      
%   end
%   
%     function selection(src,event)
%         val = c.Value;
%         disp(['Selection: ' dropboxFolderNameCell{val}]);
%     end
end
