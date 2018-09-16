function coatChannel = getCoatChannel(ExperimentType, Channel1, Channel2)

  % (MT, 2018-02-11) Added support for lattice imaging, maybe temporary -
  % FIX LATER
  coatChannel = [];

  if strcmpi(ExperimentType, 'inputoutput') || strcmpi(ExperimentType, 'lattice')

    if contains(Channel2, 'mcp', 'IgnoreCase', true) || ...
      contains(Channel2, 'pcp', 'IgnoreCase', true)
      coatChannel = 2;
    elseif contains(Channel1, 'mcp', 'IgnoreCase', true) || ...
      contains(Channel1, 'pcp', 'IgnoreCase', true)
      coatChannel = 1;
    else 
      error('No MCP or PCP channel detected. Check MovieDatabase.XLSX')
    end 

  end 

end 