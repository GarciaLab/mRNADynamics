function spotChannels = getCoatChannel(Channel1, Channel2, varargin)
    %Finds the channel that contains the stem loop coat protein (aka the
    %channel with the transcription spots
    spotChannels = [];
    
    if ~isempty(varargin)
        Channel3 = varargin{1};
    end
    
    channels = {Channel1, Channel2, Channel3};
    
    spotChannels = find(contains(channels,'CP', 'IgnoreCase',true)...
        | contains(channels,'Spot', 'IgnoreCase',true));

end 
