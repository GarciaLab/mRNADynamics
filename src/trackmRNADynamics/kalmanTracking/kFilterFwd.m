function [pdTrack,ctTrack,pdTrackSE,kalmanFilter ] = kFilterFwd(posData,InitError,MotionNoise,MeasurementNoise,timePoints)


  % initialize particle filter
  kalmanFilter =  configureKalmanFilter('ConstantAcceleration', ...
            posData(1,:), InitError, MotionNoise,MeasurementNoise);


  trackLen = size(posData,1);
  nIter = max(timePoints);
  
  % make first prediction
  pdTrack = [posData(1,:)];
  ctTrack = [posData(1,:)];
  pdTrackSE = [0 0];

  for k = 2:nIter 
    
    % make prediction
    [PredictedMeasurement,PredictedState,PredictedStateCovariance] = predict(kalmanFilter);
    % record
    pdTrack(k,:) = PredictedMeasurement';
    pdTrackSE(k,:) = [PredictedStateCovariance(1,1) PredictedStateCovariance(4,4)];  
    
    if k <=trackLen
      % correct
      ctTrack(k,:) = correct(kalmanFilter,posData(k,:)); 
    else
      ctTrack(k,:) = NaN;
    end
            
  end
%   logLTrack = sum(logLtrack);