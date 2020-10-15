function KFTrack =...
                kalmanFilterFwd(posData,InitError,MotionNoise,MeasurementNoise) 
  
  KFTrack.kfType = 'ConstantVelocity';
  if strcmp(KFTrack.kfType,'ConstantVelocity')
    nt = 2;
  elseif strcmp(KFTrack.kfType,'ConstantAcceleration')
    nt = 3;
  end
  % initialize particle filter
  KFTrack.kalmanFilter =  configureKalmanFilter(KFTrack.kfType, ...
            posData(1,:), repelem(InitError,nt), repelem(MotionNoise,nt), MeasurementNoise);               
           
  % get size 
  nDims = size(posData,2);
  nSteps = size(posData,1);
  
  % initialize arrays  
  KFTrack.priorState = NaN(nSteps,nt*nDims);
  KFTrack.priorCovariance = NaN(nt*nDims,nt*nDims,nSteps);
    
  KFTrack.posteriorState = NaN(nSteps,nt*nDims);  
  KFTrack.posteriorCovariance = NaN(nt*nDims,nt*nDims,nSteps);
  
  % set posterior values to initial values  
  KFTrack.posteriorState(1,:) = [posData(1,1) repelem(0,nt-1) posData(1,2) repelem(0,nt-1)]; 
  KFTrack.posteriorCovariance(:,:,1) = diag(repelem(InitError,nt*nDims));  
  
  for k = 2:nSteps    
    % predict
    [~,KFTrack.priorState(k,:),KFTrack.priorCovariance(:,:,k)] = predict(KFTrack.kalmanFilter);
    
    % correct
    if all(~isnan(posData(k,:)))      
       [~,KFTrack.posteriorState(k,:),KFTrack.posteriorCovariance(:,:,k)] = correct(KFTrack.kalmanFilter,posData(k,:));               
    else
      KFTrack.posteriorCovariance(:,:,k) = KFTrack.priorCovariance(:,:,k);
      KFTrack.posteriorState(k,:) = KFTrack.priorState(k,:);         
    end            
  end
  
  dimVec = 1:nt:nDims*nt;
  
  % add useful fields
  KFTrack.posData = posData;
    
  KFTrack.priorTrack = KFTrack.priorState(:,dimVec);
  KFTrack.posteriorTrack = KFTrack.posteriorState(:,dimVec);
  for d = 1:length(dimVec)
      KFTrack.priorTrackSE(:,d) = sqrt(reshape(KFTrack.priorCovariance(dimVec(d),dimVec(d),:),[],1));
      KFTrack.posteriorTrackSE(:,d) = sqrt(reshape(KFTrack.posteriorCovariance(dimVec(d),dimVec(d),:),[],1));
  end
                        