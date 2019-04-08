function coatChannel = getCoatChannel(ExperimentType, Channel1, Channel2)

  % (MT, 2018-02-11) Added support for lattice imaging, maybe temporary -
  % FIX LATER
  % (MT, 2019-04-03) Hacky fix to run input2spot data through as 2spot. 
  % Does this even need to take the ExperimentType into account?
  coatChannel = [];

  if strcmpi(ExperimentType, 'inputoutput') || strcmpi(ExperimentType, 'lattice') || strcmpi(ExperimentType, '2spot')

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
