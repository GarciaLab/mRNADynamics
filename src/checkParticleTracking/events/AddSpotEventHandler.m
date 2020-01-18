function [textInputHandler, keyInputHandler] =...
    AddSpotEventHandler(cptState, smart_add_spot, PreProcPath, FilePrefix, Prefix, robot, fake_event)
 
    function textInput(add_spot, event)
        Overlay = ancestor(add_spot, 'figure');
        figure(Overlay);

        smart_add = '{';
        
        if smart_add_spot.Value
            smart_add = '[';
        end

        PathPart1 = [PreProcPath, filesep, FilePrefix(1:end - 1), filesep, FilePrefix];
        PathPart2 = [cptState.nameSuffix, '.tif'];
        Path3 = [PreProcPath, filesep, Prefix, filesep, Prefix];
        cptState.no_clicking = true;

        [xSize, ySize, pixelSize, zStep, snippet_size, LinesPerFrame, PixelsPerLine,...
            numFrames, NDigits] = getFrameInfoParams(FrameInfo);

        [cptState.SpotFilter, cptState.Particles, cptState.Spots, cptState.PreviousParticle, cptState.CurrentParticle] = ...
            addSpot(cptState.ZoomMode, cptState.GlobalZoomMode, cptState.Particles, cptState.CurrentChannel, ...
            cptState.CurrentParticle, cptState.CurrentFrame, cptState.CurrentZ, Overlay, snippet_size, PixelsPerLine, ...
            LinesPerFrame, cptState.Spots, cptState.ZSlices, PathPart1, PathPart2, Path3, cptState.FrameInfo, pixelSize, ...
            cptState.SpotFilter,smart_add, xSize, ySize, NDigits, ...
            cptState.coatChannel, cptState.UseHistoneOverlay, cptState.schnitzcells, cptState.nWorkers, cptState.plot3DGauss);
        
        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);
        cptState.no_clicking = false;
    end

    function keyInput(cc)

    end

    textInputHandler = @textInput;
    keyInputHandler = @keyInput;
end
