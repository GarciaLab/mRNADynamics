function trackingInfo = parseTrackingOptions(liveExperiment, Spots, useHistone)

  FrameInfo = getFrameInfo(liveExperiment);
  Tres = median(diff([FrameInfo.Time]));
  
  % get number of channels
  trackingInfo.NCh = length(Spots);
  
  % set type of kalman filter
  trackingInfo.kfType = 'ConstantVelocity';
  
  % initialize linking hyperparameters
  trackingInfo.matchCostMax = repelem(realmax,trackingInfo.NCh);  
  trackingInfo.maxUnobservedFrames = repelem(Inf,trackingInfo.NCh);  
  trackingInfo.spotsPerNucleus = repelem(Inf,trackingInfo.NCh);  
  
  if ~useHistone && ~ismember(liveExperiment.experimentType,{'inputoutput'})
      trackingInfo.matchCostMax = repelem(-2*log(1e-7),trackingInfo.NCh);  
      trackingInfo.maxUnobservedFrames = repelem(Inf,trackingInfo.NCh);
      
  elseif ~useHistone && ismember(liveExperiment.experimentType,{'inputoutput'})
      % Figure out which channel is the cluster channel
      spotsChannels = liveExperiment.spotChannels;
      inputChannels = liveExperiment.inputChannels;  
      ClusterFilter = ismember(spotsChannels,inputChannels);  

      % Place reasonable limit on link cost for clusters
      trackingInfo.matchCostMax = repelem(-2*log(1e-7),trackingInfo.NCh); % this seems to work reasonably well 

      % Place ceiling on number of frames a cluster an go undetected
      % before being "capped"
      trackingInfo.maxUnobservedFrames(ClusterFilter) = ceil(3*60 / Tres);
      
  elseif ismember(liveExperiment.experimentType,{'1spot'}) 
      trackingInfo.spotsPerNucleus = 1;    
      
  elseif ismember(liveExperiment.experimentType,{'2spot'}) 
      trackingInfo.spotsPerNucleus = 2;    
      
  elseif ismember(liveExperiment.experimentType,{'2spot2color'}) && trackingInfo.NCh==2
      trackingInfo.spotsPerNucleus = [1,1];      
      
  elseif ismember(liveExperiment.experimentType,{'inputoutput'}) && trackingInfo.NCh <= 3
    
      % Figure out which channel is the cluster channel
      spotsChannels = liveExperiment.spotChannels;
      inputChannels = liveExperiment.inputChannels;  
      ClusterFilter = ismember(spotsChannels,inputChannels);

      % Cluster channel not limited in spotsPerNucleus        
      trackingInfo.spotsPerNucleus(~ClusterFilter) = 1;

      % Place reasonable limit on link cost for clusters
      trackingInfo.matchCostMax(ClusterFilter) = -2*log(1e-7); % this seems to work reasonably well 

      % Place ceiling on number of frames a cluster an go undetected
      % before being "capped"
      trackingInfo.maxUnobservedFrames(ClusterFilter) = ceil(3*60 / Tres); 
      
  else
      error(['''',liveExperiment.experimentType,''' liveExperiment.experimentType not supported by track04StitchTracks'])
  end