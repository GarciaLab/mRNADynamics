function KFTrack = kalmanFilterFwd(posData,kalmanOptions) 
  
  KFTrack.type = kalmanOptions.type;
  if strcmp(KFTrack.type,'ConstantVelocity')
    nt = 2;
  elseif strcmp(KFTrack.type,'ConstantAcceleration')
    nt = 3;
  end
  
  % initialize particle filter 
  KFTrack.kalmanFilter =  configureKalmanFilter(kalmanOptions.type, ...
            posData(1,:), kalmanOptions.InitialError, kalmanOptions.MotionNoise, kalmanOptions.MeasurementNoise);               
           
  % get size 
  nDims = size(posData,2);
  nSteps = size(posData,1);
  
  % initialize arrays  
  KFTrack.logL = NaN(nSteps,1);
  
  KFTrack.priorState = NaN(nSteps,nt*nDims);
  KFTrack.priorCovariance = NaN(nt*nDims,nt*nDims,nSteps);
    
  KFTrack.posteriorState = NaN(nSteps,nt*nDims);  
  KFTrack.posteriorCovariance = NaN(nt*nDims,nt*nDims,nSteps);
  
  % set posterior values to initial values  
  for i = 1:nDims
      ind = (i-1)*nt+1;
      KFTrack.posteriorState(1,ind:ind+nt-1) = [posData(1,i) repelem(0,nt-1)]; 
  end
  KFTrack.posteriorCovariance(:,:,1) = diag(repelem(kalmanOptions.InitialError,nDims));  
  
  for k = 2:nSteps    
    % predict
    [~,KFTrack.priorState(k,:),KFTrack.priorCovariance(:,:,k)] = predict(KFTrack.kalmanFilter);
    
    % update
    if all(~isnan(posData(k,:))) 
       KFTrack.logL(k) = -distance(KFTrack.kalmanFilter,posData(k,:));
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
                        