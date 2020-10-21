function KFTrack = kalmanFilterBkd(KFTrack)
      

  KFTrack.kfType = 'ConstantVelocity';
  if strcmp(KFTrack.kfType,'ConstantVelocity')
    nt = 2;
  elseif strcmp(KFTrack.kfType,'ConstantAcceleration')
    nt = 3;
  end
  
  % Initialize
  KFTrack.smoothedState = NaN(size(KFTrack.priorState));
  KFTrack.smoothedCovariance = NaN(size(KFTrack.priorCovariance));

  last_i = size(KFTrack.priorTrack,1);
  F = KFTrack.kalmanFilter.StateTransitionModel;

  KFTrack.smoothedState(last_i,:) = KFTrack.posteriorState(last_i,:);
  KFTrack.smoothedCovariance(:,:,last_i) = KFTrack.posteriorCovariance(:,:,last_i);

  for it = last_i-1:-1:1
      C = KFTrack.posteriorCovariance(:,:,it) * F' / KFTrack.priorCovariance(:,:,it+1);
      KFTrack.smoothedState(it,:) = KFTrack.posteriorState(it,:)'   + C * (KFTrack.smoothedState(it+1,:)'   - KFTrack.priorState(it+1,:)');
      KFTrack.smoothedCovariance(:,:,it) = KFTrack.posteriorCovariance(:,:,it) + C * (...
        KFTrack.smoothedCovariance(:,:,it+1) - KFTrack.priorCovariance(:,:,it+1)) * C';
  end
  
  dimVec = 1:nt:size(KFTrack.smoothedState,2);
  % add useful fields      
  KFTrack.smoothedTrack = KFTrack.smoothedState(:,dimVec);  
  for d = 1:length(dimVec)
      KFTrack.smoothedTrackSE(:,d) = sqrt(reshape(KFTrack.smoothedCovariance(dimVec(d),dimVec(d),:),[],1));      
  end
      
      
      