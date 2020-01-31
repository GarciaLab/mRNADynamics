function keyInputHandler = GeneralEventHandler(cptState, DataFolder, DropboxFolder, FilePrefix, NChannels)
 
    function keyInput(cc)
        if cc == 'f'
            [cptState.Particles, cptState.schnitzcells] = redoTracking(DataFolder, ...
                cptState.UseHistoneOverlay, cptState.FrameInfo, DropboxFolder, FilePrefix, cptState.schnitzcells, ...
                cptState.Particles, NChannels, cptState.CurrentChannel, cptState.numParticles());
        elseif cc == 's'
            saveChanges(NChannels, cptState.Particles, cptState.Spots, cptState.SpotFilter, DataFolder, ...
                cptState.FrameInfo, cptState.UseHistoneOverlay, FilePrefix, ...
                cptState.schnitzcells, DropboxFolder);
        elseif cc == '~'
            % Switch projection mode
            cptState.projectionMode = chooseProjection;
            disp(['projectionMode : ' cptState.projectionMode])
        end
    end

    keyInputHandler = @keyInput;
end
