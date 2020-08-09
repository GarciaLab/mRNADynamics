function keyInputHandler = AddSpotEventHandler(cptState, Prefix)
%Add particle and all of its shadows to cptState.Spots.

    function keyInput(cc)
        if cc == '[' | cc == '{' %#ok<*OR2>
            liveExperiment = LiveExperiment(Prefix);   
            movieMat = getMovieMat(liveExperiment);
            imStack = double(movieMat(:, :, :, cptState.CurrentFrame, cptState.CurrentChannel));
            
            [cptState.ParticleStitchInfo, cptState.SpotFilter, cptState.Particles, cptState.Spots, cptState.PreviousParticle,...
                cptState.CurrentParticle, cptState.ZoomMode, cptState.GlobalZoomMode] = ...
                ...
                addSpot(cptState.ZoomMode, cptState.GlobalZoomMode, cptState.Particles, cptState.ParticleStitchInfo, cptState.CurrentChannelIndex, ...
                ...
                cptState.CurrentParticle, cptState.CurrentFrame, cptState.CurrentZ, ...
                cptState.Spots, cptState.SpotFilter, cc, ...
                Prefix, cptState.UseHistoneOverlay, cptState.schnitzcells, cptState.nWorkers, cptState.plot3DGauss, imStack);
              
            % write file instructing pipeline to re-run tracking
            rerunParticleTracking = struct;
            save([liveExperiment.resultsFolder, 'rerunParticleTracking.mat'],'rerunParticleTracking');
        end
    end

keyInputHandler = @keyInput;
end
