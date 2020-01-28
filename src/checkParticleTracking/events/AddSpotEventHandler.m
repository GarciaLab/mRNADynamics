function keyInputHandler = AddSpotEventHandler(cptState, PreProcPath, ProcPath, Prefix)
%Add particle and all of its shadows to cptState.Spots.

    function keyInput(cc)
        if cc == '[' | cc == '{' %#ok<*OR2>
            PathPart1 = [PreProcPath, filesep, Prefix, filesep, Prefix, '_'];
            PathPart2 = [cptState.nameSuffix, '.tif'];
            Path3 = [PreProcPath, filesep, Prefix, filesep, Prefix];
            
            [xSize, ySize, pixelSize, zStep, snippet_size, LinesPerFrame, PixelsPerLine,...
                nFrames, NSlices] = getFrameInfoParams(cptState.FrameInfo);
            
            %See how  many frames we have and adjust the index size of the files to load accordingly
            if nFrames < 1E3
                NDigits = 3;
            elseif nFrames < 1E4
                NDigits = 4;
            else
                error('No more than 10,000 frames supported.')
            end
            
            [cptState.SpotFilter, cptState.Particles, cptState.Spots, cptState.PreviousParticle, cptState.CurrentParticle, cptState.ZoomMode, cptState.GlobalZoomMode] = ...
                addSpot(cptState.ZoomMode, cptState.GlobalZoomMode, cptState.Particles, cptState.CurrentChannel, ...
                cptState.CurrentParticle, cptState.CurrentFrame, cptState.CurrentZ, snippet_size, PixelsPerLine, ...
                LinesPerFrame, cptState.Spots, cptState.ZSlices, PathPart1, PathPart2, Path3, cptState.FrameInfo, pixelSize, ...
                cptState.SpotFilter, cc, xSize, ySize, NDigits, ...
                Prefix, PreProcPath, ProcPath, cptState.coatChannel, cptState.UseHistoneOverlay, cptState.schnitzcells, cptState.nWorkers, cptState.plot3DGauss);
        end
    end

keyInputHandler = @keyInput;
end
