function [textInputHandler, keyInputHandler] =...
    AddSpotEventHandler(cptState, smart_add_spot, PreProcPath, ProcPath, Prefix, robot, fake_event)
 %Add particle and all of its shadows to cptState.Spots.
 
    function doAddSpot(PreProcPath, ProcPath,Prefix, cc)
        PathPart1 = [PreProcPath, filesep, Prefix, filesep, Prefix, '_'];
        PathPart2 = [cptState.nameSuffix, '.tif'];
        Path3 = [PreProcPath, filesep, Prefix, filesep, Prefix];
        
        [xSize, ySize, pixelSize, zStep, snippet_size, LinesPerFrame, PixelsPerLine,...
        numFrames, NDigits] = getFrameInfoParams(cptState.FrameInfo);

        [cptState.SpotFilter, cptState.Particles, cptState.Spots, cptState.PreviousParticle, cptState.CurrentParticle, cptState.ZoomMode, cptState.GlobalZoomMode] = ...
            addSpot(cptState.ZoomMode, cptState.GlobalZoomMode, cptState.Particles, cptState.CurrentChannel, ...
            cptState.CurrentParticle, cptState.CurrentFrame, cptState.CurrentZ, snippet_size, PixelsPerLine, ...
            LinesPerFrame, cptState.Spots, cptState.ZSlices, PathPart1, PathPart2, Path3, cptState.FrameInfo, pixelSize, ...
            cptState.SpotFilter, cc, xSize, ySize, NDigits, ...
            Prefix, PreProcPath, ProcPath, cptState.coatChannel, cptState.UseHistoneOverlay, cptState.schnitzcells, cptState.nWorkers, cptState.plot3DGauss);
    end

    function textInput(add_spot, event)
        Overlay = ancestor(add_spot, 'figure');
        figure(Overlay);

        cc = '{';
        
        if smart_add_spot.Value
            cc = '[';
        end

        cptState.no_clicking = true;

        doAddSpot(PreProcPath, ProcPath, Prefix, cc);
        
        robot.keyPress(fake_event);
        robot.keyRelease(fake_event);
        cptState.no_clicking = false;
    end

    function keyInput(cc)
        if cc == '[' | cc == '{' %#ok<*OR2>
            doAddSpot(PreProcPath, ProcPath, Prefix, cc);
        end
    end

    textInputHandler = @textInput;
    keyInputHandler = @keyInput;
end
