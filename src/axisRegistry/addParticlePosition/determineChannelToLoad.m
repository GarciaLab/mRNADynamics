function ChannelToLoad = determineChannelToLoad(SelectChannel, Channels)

    if SelectChannel
        list = string(Channels);
        [indx,tf] = listdlg('PromptString','Select the channel to use for alignment:','ListString',list);
        ChannelToLoad = indx;
    else
        % From now, we will use a better way to define the channel for
        % alignment (used for cross-correlation).
        % Find channels with ":Nuclear"
        ChannelToLoadTemp= contains([Channels(1),Channels(2),Channels(3)],'nuclear','IgnoreCase',true);
        
        % Define the Channel to load, for calculating the cross-correlation
        % In future, we can think about combining multiple channels for
        % calculating the cross-correlation to get more accurate estimates.
        % For now, let's pick only one channel for this. For multiple
        % channels, let's first pick the first channel. This can be fine in
        % most cases, since we normally use lower wavelength for sth we
        % care more, or we get better signal from those.
        if sum(ChannelToLoadTemp) && sum(ChannelToLoadTemp)==1
            ChannelToLoad=find(ChannelToLoadTemp);
        elseif sum(ChannelToLoadTemp) && length(ChannelToLoadTemp)>=2
            ChannelToLoad=find(ChannelToLoadTemp);
            ChannelToLoad = ChannelToLoad(1);
        else
            error('No histone channel found. Was it defined in MovieDatabase as :Nuclear or :InvertedNuclear?')
        end
        
    end
end
