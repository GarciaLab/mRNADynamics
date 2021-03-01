function keyInputHandler = GeneralEventHandler(cptState, DataFolder, DropboxFolder, FilePrefix, NChannels)
 
    function keyInput(cc)
        if cc == 'f'
            [cptState.Particles, cptState.schnitzcells] = redoTracking(DataFolder, ...
                cptState.UseHistoneOverlay, cptState.FrameInfo, DropboxFolder, FilePrefix, cptState.schnitzcells, ...
                cptState.Particles, NChannels, cptState.CurrentChannelIndex, cptState.numParticles());
        elseif cc == 's'
            %disp('Temporarily disabled save.')
            saveChanges(NChannels, cptState, DataFolder, FilePrefix, DropboxFolder);
        elseif cc == '~'
            % Switch projection mode
            cptState.projectionMode = chooseProjection;
            disp(['projectionMode : ' cptState.projectionMode])
        
        elseif cc == '!'
            % Increase contrast in the Overlay figure
            if isempty(cptState.DisplayRangeSpot)
                cptState.DisplayRangeSpot = [min(min(cptState.ImageMat)), max(max(cptState.ImageMat)) / 1.5];
            else
                cptState.DisplayRangeSpot = [cptState.DisplayRangeSpot(1), cptState.DisplayRangeSpot(2) / 1.5];
            end
            
            disp('increased spot contrast');
        
        elseif cc == '@'
            % Decrease spot channel contrast
            cptState.DisplayRangeSpot = [min(cptState.ImageMat(:)), max(cptState.ImageMat(:)) * 1.5];
            
            disp('decreased spot contrast');
            
        elseif cc == '0' %Debugging mode
            keyboard;
        end
    end

    keyInputHandler = @keyInput;
end
