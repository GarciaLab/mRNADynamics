%%
% Given a cell array in the form of key-value pairs,
% returns the value corresponding to the property "propertyName".
%% 
function value = getConfigValue(configuration, propertyName)
  configColumn = configuration(:, 1);
  configIndexArray = strfind(configColumn, propertyName);
  configIndex = find(not(cellfun('isempty', configIndexArray)));

  if isempty(configIndex)
    error(['Property ''', propertyName, ''' not found in configuration. ',...
      'Check ComputerFolders and/or MovieDatabase.'])
  end

  valueCell = configuration(configIndex, 2);
  value = valueCell{1};
end
