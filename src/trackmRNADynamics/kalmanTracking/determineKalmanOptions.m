function kalmanOptions = determineKalmanOptions(liveExperiment,trackingInfo, Spots)
    
    if strcmp(trackingInfo.kfType,'ConstantVelocity')
      nDims = 2;
    elseif strcmp(trackingInfo.kfType,'ConstantAcceleration')
      nDims = 3;
    end
    kalmanOptions.type = trackingInfo.kfType;
    
    %% Set noise parameters
    kalmanOptions.MeasurementNoise = 0.1/liveExperiment.pixelSize_um; 
    kalmanOptions.MotionNoise = repelem(kalmanOptions.MeasurementNoise,nDims);
    kalmanOptions.InitialError = repelem(kalmanOptions.MeasurementNoise,nDims);
    
    
    %% Estimate adjustment factor for fluorescence field
    if ~trackingInfo.use3DInfo
        fluoVec = [];
        for s = 1:length(Spots)
            Fits = Spots(s).Fits;
            for f = 1:length(Fits)
                fluoVec = [fluoVec Fits(f).FixedAreaIntensity3/3 - Fits(f).FixedAreaIntensity(Fits(f).z==Fits(f).brightestZ)];
            end
        end
        kalmanOptions.fluoFactor = sqrt(mean(fluoVec.^2))*5; % factor of 5 functions
        
    else
        fluoVec = [];
        for s = 1:length(Spots)
            Fits = Spots(s).Fits;
            if ~isempty(Fits)
                fluoVec = [fluoVec [Fits.gauss3DIntensityCI95]];
            end
        end
        kalmanOptions.fluoFactor = sqrt(mean(fluoVec.^2))*5;
    end
    
    %% Should we include distance from nucleus center?
    kalmanOptions.useNuclei = trackingInfo.useHistone;
    
    if kalmanOptions.useNuclei
        kalmanOptions.measurementFields = {'xPos', 'yPos', 'zPos', 'Fluo', 'nucleusDist'};
    else
        kalmanOptions.measurementFields = {'xPos', 'yPos', 'zPos', 'Fluo'};
    end