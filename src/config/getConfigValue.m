%%
% Given a cell array in the form of key-value pairs,
% returns the value corresponding to the property "propertyName".
%% 
function value = getConfigValue(configuration, propertyName)

arguments
    configuration string
    propertyName string
end

configColumn = configuration(:, 1);
configIndex = find(configColumn == propertyName); 

  if isempty(configIndex)
    error(['Property ''', char(propertyName), ''' not found in configuration. ',...
      'Check ComputerFolders and/or MovieDatabase.'])
  end
  
  %cast to char for compatibility
  value = char(strtrim(configuration(configIndex, 2))); 
end
