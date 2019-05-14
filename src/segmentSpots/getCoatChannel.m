function spotChannels = getCoatChannel(Channel1, Channel2, varargin)
    %Finds the channel that contains the stem loop coat protein (aka the
    %channel with the transcription spots
    spotChannels = [];
    
    if ~isempty(varargin)
        Channel3 = varargin{1};
    end
    
    spotChannels = find(contains([Channel1,Channel2,Channel3],'CP') || contains([Channel1,Channel2,Channel3],'Spot'));

end 
