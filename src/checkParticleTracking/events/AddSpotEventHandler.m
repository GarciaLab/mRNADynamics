function keyInputHandler = AddSpotEventHandler(cptState, PreProcPath, ProcPath, Prefix, movieMat)
%Add particle and all of its shadows to cptState.Spots.

    function keyInput(cc)
        if cc == '[' | cc == '{' %#ok<*OR2>
            PathPart1 = [PreProcPath, filesep, Prefix, filesep, Prefix, '_'];
            PathPart2 = [cptState.nameSuffix, '.tif'];
            Path3 = [PreProcPath, filesep, Prefix, filesep, Prefix];
            
            [xSize, ySize, pixelSize, zStep, snippet_size,...
                nFrames, NSlices, NDigits] = getFrameInfoParams(cptState.FrameInfo);
          
           
            
            [cptState.SpotFilter, cptState.Particles, cptState.Spots, cptState.PreviousParticle, cptState.CurrentParticle, cptState.ZoomMode, cptState.GlobalZoomMode] = ...
                addSpot(cptState.ZoomMode, cptState.GlobalZoomMode, cptState.Particles, cptState.CurrentChannel, ...
                cptState.CurrentParticle, cptState.CurrentFrame, cptState.CurrentZ, snippet_size, xSize, ...
                ySize, cptState.Spots, cptState.ZSlices, PathPart1, PathPart2, Path3, cptState.FrameInfo, pixelSize, ...
                cptState.SpotFilter, cc, xSize, ySize, NDigits, ...
                Prefix, PreProcPath, ProcPath, cptState.coatChannel, cptState.UseHistoneOverlay, cptState.schnitzcells, cptState.nWorkers, cptState.plot3DGauss, movieMat);
        end
    end

keyInputHandler = @keyInput;
end
