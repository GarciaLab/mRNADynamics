function keyInputHandler = AddSpotEventHandler(cptState, Prefix)
%Add particle and all of its shadows to cptState.Spots.

    function keyInput(cc)
        if cc == '[' | cc == '{' | cc == '}' %#ok<*OR2>
            if cc == '}'
                cptState = removeCurrentFrame(cptState);
            end
            movieMat = getMovieMat(LiveExperiment(Prefix));
            if ~isempty(movieMat)
                imStack = double(movieMat(:, :, :, cptState.CurrentFrame, cptState.CurrentChannel));
            else
                imStack = getMovieFrame(LiveExperiment(Prefix), cptState.CurrentFrame, cptState.CurrentChannel);
            end
            
            [cptState.SpotFilter, cptState.Particles, cptState.Spots, cptState.PreviousParticle,...
                cptState.CurrentParticle, cptState.ZoomMode, cptState.GlobalZoomMode] = ...
                ...
                addSpot(cptState.ZoomMode, cptState.GlobalZoomMode, cptState.Particles, cptState.CurrentChannelIndex, ...
                ...
                cptState.CurrentParticle, cptState.CurrentFrame, cptState.CurrentZ, ...
                cptState.Spots, cptState.SpotFilter, cc, ...
                Prefix, cptState.UseHistoneOverlay, cptState.schnitzcells, cptState.nWorkers, cptState.plot3DGauss, imStack);
        end
    end

keyInputHandler = @keyInput;
end
