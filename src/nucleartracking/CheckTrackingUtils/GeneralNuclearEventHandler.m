function keyInputHandler = GeneralNuclearEventHandler(cntState, DataFolder, DropboxFolder, FilePrefix)
 
    function keyInput(cc)

        if cc == 's'
            saveNuclearChanges(cntState, DataFolder, FilePrefix, DropboxFolder);
        elseif cc == '~'
            % Switch projection mode
            cntState.projectionMode = chooseProjection;
            disp(['projectionMode : ' cntState.projectionMode])
        
        elseif cc == '!'
            % Increase contrast in the Overlay figure
            if isempty(cntState.DisplayRangeSpot)
                cntState.DisplayRangeSpot = [min(min(cntState.ImageMat)), max(max(cntState.ImageMat)) / 1.5];
            else
                cntState.DisplayRangeSpot = [cntState.DisplayRangeSpot(1), cntState.DisplayRangeSpot(2) / 1.5];
            end
            
            disp('increased spot contrast');
        
        elseif cc == '@'
            % Decrease spot channel contrast
            cntState.DisplayRangeSpot = [min(cntState.ImageMat(:)), max(cntState.ImageMat(:)) * 1.5];
            
            disp('decreased spot contrast');
            
        elseif cc == '0' %Debugging mode
            keyboard;
        end
    end

    keyInputHandler = @keyInput;
end
