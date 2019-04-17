function coatChannel = getCoatChannel(ExperimentType, Channel1, Channel2)
    %Finds the channel that contains the stem loop coat protein (aka the
    %channel with the transcription spots)
    coatChannel = [];

    if strcmpi(ExperimentType, 'input') 
         coatChannel = [];
    elseif contains(Channel2, 'mcp', 'IgnoreCase', true) || ...
      contains(Channel2, 'pcp', 'IgnoreCase', true)
      coatChannel = 2;
    elseif contains(Channel1, 'mcp', 'IgnoreCase', true) || ...
      contains(Channel1, 'pcp', 'IgnoreCase', true)
      coatChannel = 1; 
    else
      error('No MCP or PCP channel detected. Check MovieDatabase.csv')
    end 

  end 

end 
