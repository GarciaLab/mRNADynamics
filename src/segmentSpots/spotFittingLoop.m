function SpotsFr = spotFittingLoop(FitsFr, liveExperiment, imStack, spotDims)
    SpotsFr.Fits = FitsFr;
    for i = 1:length(SpotsFr.Fits)
        SpotsFr = fitSnip3D(SpotsFr, i, liveExperiment, imStack, spotDims);    
    end    