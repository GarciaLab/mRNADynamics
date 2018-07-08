%Figure out the different channels
function [coatChannel, histoneChannel, fiducialChannel] = obtainZeissChannels(Channel1, Channel2)
  if ~isempty(strfind(Channel1{1}, 'MCP'))
    coatChannel = 1;
  elseif strfind(Channel1{1}, 'His')
    histoneChannel = 1;
  else
    error('LSM Mode error: Channel name not recognized. Check MovieDatabase')
  end

  if ~isempty(strfind(Channel2{1}, 'MCP'))
    coatChannel = 2;
  elseif strfind(Channel2{1}, 'His')
    histoneChannel = 2;
  else
    error('LSM Mode error: Channel name not recognized. Check MovieDatabase')
  end

  fiducialChannel = histoneChannel;
end
