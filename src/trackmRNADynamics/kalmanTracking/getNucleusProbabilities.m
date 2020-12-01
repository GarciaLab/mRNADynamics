function Particles = getNucleusProbabilities(liveExperiment,Particles,FrameInfo,trackingOptions)

    %%% load nucleus probability maps if they exist
    nucleusProbDirRaw = [liveExperiment.procFolder 'nucleusProbabilityMaps' filesep];
    nucleusProbDirFinal = [liveExperiment.procFolder 'nucleusProbabilityMapsFull' filesep];
    probFiles = dir([nucleusProbDirFinal '*.tif']);
    probFilesRaw = dir([nucleusProbDirRaw '*.tif']);
    hasFullProbFiles = length(probFiles) == length(FrameInfo);

    if hasFullProbFiles
        disp('Incorporating nucleus probabilities...')
        Particles = addNucleusProbabilities(liveExperiment, trackingOptions, FrameInfo, Particles);      
    else
        if exist([nucleusProbDirRaw 'nucleusInfo.mat'], 'file')
            load([nucleusProbDirRaw 'nucleusInfo.mat'], 'nucleusInfo')
            if length(nucleusInfo.originalFileNames) == length(probFilesRaw)
                compileNuclearTiffStacks(liveExperiment.Prefix);
                disp('Incorporating nucleus probabilities...')
                Particles = addNucleusProbabilities(liveExperiment, trackingOptions, FrameInfo, Particles);      
            else
                warning('No nucleus probabilities found. Nucleus boundaries will not be used to assess particle quality.')
                for Channel = 1:trackingOptions.NCh
                  for  p = 1:length(Particles{Channel})
                      Particles{Channel}(p).nucleusProbability = ones(size(Particles{Channel}(p).Frame));
                  end
                end
            end
        else
            warning('No nucleus probabilities found. Nucleus boundaries will not be used to assess particle quality.')
            for Channel = 1:trackingOptions.NCh
                for  p = 1:length(Particles{Channel})
                    Particles{Channel}(p).nucleusProbability = ones(size(Particles{Channel}(p).Frame));
                end
            end
        end
    end