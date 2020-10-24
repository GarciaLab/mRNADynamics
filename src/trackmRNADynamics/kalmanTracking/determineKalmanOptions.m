function kalmanOptions = determineKalmanOptions(liveExperiment,kfType)
    
    if strcmp(kfType,'ConstantVelocity')
      nDims = 2;
    elseif strcmp(kfType,'ConstantAcceleration')
      nDims = 3;
    end
    kalmanOptions.type = kfType;
    % NL: these parameters generally work pretty well
    kalmanOptions.MeasurementNoise = 0.1/liveExperiment.pixelSize_um; 
    kalmanOptions.MotionNoise = repelem(kalmanOptions.MeasurementNoise,nDims);%e-3; % NL: this seems to work well
    kalmanOptions.InitialError = repelem(kalmanOptions.MeasurementNoise,nDims);