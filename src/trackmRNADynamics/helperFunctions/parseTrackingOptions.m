function trackingInfo = parseTrackingOptions(liveExperiment, Spots, useHistone)

  FrameInfo = getFrameInfo(liveExperiment);
  Tres = median(diff([FrameInfo.Time]));
  
  % should we use 3D spot info? NL: leaving this 2D until 3D fits are
  % improved
  trackingInfo.use3DInfo = false;
  
  % designate which observables to use for tracking (stick with position
  % for now)
  trackingInfo.trackingIndices = 1:3;
  
  % record nucleus flag
  trackingInfo.useHistone = useHistone;
  
  % get number of channels
  trackingInfo.NCh = length(Spots);
  
  % number of frames 
  trackingInfo.nFrames = length(FrameInfo);
  
  % set number of frames to extrapolate before and after last detections
  trackingInfo.nExtrapFrames = ceil(2*60/Tres);
  
  % set type of kalman filter
  trackingInfo.kfType = 'ConstantVelocity';
  
  % initialize linking hyperparameters
  trackingInfo.matchCostMax = repelem(realmax,trackingInfo.NCh);  
  trackingInfo.maxUnobservedFrames = repelem(Inf,trackingInfo.NCh);  
  trackingInfo.spotsPerNucleus = repelem(Inf,trackingInfo.NCh);  
  trackingInfo.targetMatchFrac = repelem(0.99,trackingInfo.NCh); 
  
  if ~useHistone && ~ismember(liveExperiment.experimentType,{'inputoutput'})
      trackingInfo.matchCostMax = repelem(-3*log(3e-5),trackingInfo.NCh);  % NL: this works pretty well. Need better way to estimate this dynamically
      trackingInfo.maxUnobservedFrames = repelem(Inf,trackingInfo.NCh);
      
  elseif ~useHistone && ismember(liveExperiment.experimentType,{'inputoutput'})
      % Figure out which channel is the cluster channel
      spotsChannels = liveExperiment.spotChannels;
      inputChannels = liveExperiment.inputChannels;  
      ClusterFilter = ismember(spotsChannels,inputChannels);  

      % Place reasonable limit on link cost for clusters
      trackingInfo.matchCostMax = repelem(-3*log(3e-5),trackingInfo.NCh); % this seems to work reasonably well 

      % Place ceiling on number of frames a cluster an go undetected
      % before being "capped"
      trackingInfo.maxUnobservedFrames(ClusterFilter) = ceil(3*60 / Tres);
      
      % adjust match fraction for clusters
      trackingInfo.targetMatchFrac(ClusterFilter) = 0.9;
      
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
      trackingInfo.matchCostMax(ClusterFilter) = -4*log(3e-5); % this seems to work reasonably well 

      % Place ceiling on number of frames a cluster an go undetected
      % before being "capped"
      trackingInfo.maxUnobservedFrames(ClusterFilter) = ceil(3*60 / Tres); 
      
      % adjust match fraction for clusters
      trackingInfo.targetMatchFrac(ClusterFilter) = 0.9;
  else
      error(['''',liveExperiment.experimentType,''' liveExperiment.experimentType not supported by track04StitchTracks'])
  end