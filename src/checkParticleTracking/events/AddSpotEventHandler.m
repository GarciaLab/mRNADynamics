function keyInputHandler = AddSpotEventHandler(cptState, Prefix)
%Add particle and all of its shadows to cptState.Spots.

    function keyInput(cc)
        if cc == '[' | cc == '{' %#ok<*OR2>
            liveExperiment = LiveExperiment(Prefix);   
            movieMat = getMovieMat(liveExperiment);
            imStack = double(movieMat(:, :, :, cptState.CurrentFrame, cptState.CurrentChannel));
          
            cptState = addSpot(cptState, cc, Prefix, imStack);
              
            % write file instructing pipeline to re-run tracking
            rerunParticleTracking = struct;
            save([liveExperiment.resultsFolder, 'rerunParticleTracking.mat'],'rerunParticleTracking');
        end
    end

keyInputHandler = @keyInput;
end
