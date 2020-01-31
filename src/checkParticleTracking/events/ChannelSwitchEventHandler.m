 function keyInputHandler = ChannelSwitchEventHandler(cptState, NChannels, coatChannels, UseHistoneOverlay)
    %NChannels, UseHistoneOverlay are never modified in the command loop, it is safe to pass them by-value here

    function keyInput(cc)
        if cc == '8' & NChannels > 1 %#ok<*AND2>
            [cptState.CurrentChannel, cptState.PreviousChannel, cptState.coatChannel, cptState.CurrentParticle] = ...
                switchChannels(cptState.CurrentChannel, cptState.CurrentParticle, cptState.Particles, ...
                UseHistoneOverlay, coatChannels, NChannels);
        end
    end

    keyInputHandler = @keyInput;
end

