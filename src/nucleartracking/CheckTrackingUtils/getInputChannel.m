function inputChannels = getInputChannel(Channel1, Channel2, varargin)
    %Finds the channel that contains the stem loop coat protein (aka the
    %channel with the transcription spots)
    inputChannels = [];
    
    if ~isempty(varargin)
        Channel3 = varargin{1};
    end
    
    
    if iscell(Channel1), Channel1 = Channel1{1}; end
    if iscell(Channel2), Channel2 = Channel2{1}; end
    if exist('Channel3', 'var') && iscell(Channel3), Channel3 = Channel3{1}; end


    channels = {Channel1, Channel2, Channel3};

    
    inputChannels = find(contains(channels,'input', 'IgnoreCase',true));

end 
