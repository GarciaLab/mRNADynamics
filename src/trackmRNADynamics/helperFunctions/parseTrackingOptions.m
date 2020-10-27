function trackingOptions = parseTrackingOptions(Spots, liveExperiment, trackingOptions)

  FrameInfo = getFrameInfo(liveExperiment);
  Tres = median(diff([FrameInfo.Time]));
  
  % should we use 3D spot info? NL: leaving this 2D until 3D fits are
  % improved
  trackingOptions.use3DInfo = false;
  
  % designate which observables to use for tracking (stick with position
  % for now)
  trackingOptions.trackingIndices = 1:3;  
  
  % set default linking cost
  trackingOptions.matchCostDefault = -3*log(3e-5);
  % get number of channels
  trackingOptions.NCh = length(Spots);
  
  % number of frames 
  trackingOptions.nFrames = length(FrameInfo);
  
  % set number of frames to extrapolate before and after last detections
  trackingOptions.nExtrapFrames = ceil(2*60/Tres);
  
  % set type of kalman filter
  trackingOptions.kfType = 'ConstantVelocity';
  
  % initialize linking hyperparameters
  trackingOptions.matchCostMaxForward = repelem(realmax,trackingOptions.NCh);  
  trackingOptions.matchCostMaxBackward = repelem(realmax,trackingOptions.NCh);  
  trackingOptions.maxUnobservedFrames = repelem(Inf,trackingOptions.NCh);  
  trackingOptions.spotsPerNucleus = repelem(Inf,trackingOptions.NCh);  
  trackingOptions.targetMatchFrac = repelem(0.99,trackingOptions.NCh); 
  
  if ~trackingOptions.useHistone && ~ismember(liveExperiment.experimentType,{'inputoutput'})
      trackingOptions.matchCostMax = repelem(trackingOptions.matchCostDefault,trackingOptions.NCh);  % NL: this works pretty well. Need better way to estimate this dynamically
      trackingOptions.maxUnobservedFrames = repelem(Inf,trackingOptions.NCh);
      
  elseif ~trackingOptions.useHistone && ismember(liveExperiment.experimentType,{'inputoutput'})
      % Figure out which channel is the cluster channel
      spotsChannels = liveExperiment.spotChannels;
      inputChannels = liveExperiment.inputChannels;  
      ClusterFilter = ismember(spotsChannels,inputChannels);  

      % Place reasonable limit on link cost for clusters
      trackingOptions.matchCostMax = repelem(-3*log(3e-5),trackingOptions.NCh); % this seems to work reasonably well 

      % Place ceiling on number of frames a cluster an go undetected
      % before being "capped"
      trackingOptions.maxUnobservedFrames(ClusterFilter) = ceil(3*60 / Tres);
      
      % adjust match fraction for clusters
      trackingOptions.targetMatchFrac(ClusterFilter) = 0.9;
      
  elseif ismember(liveExperiment.experimentType,{'1spot'}) 
      trackingOptions.spotsPerNucleus = 1;    
      
  elseif ismember(liveExperiment.experimentType,{'2spot'}) 
      trackingOptions.spotsPerNucleus = 2;    
      trackingOptions.matchCostMaxForward = trackingOptions.matchCostDefault;
      
  elseif ismember(liveExperiment.experimentType,{'2spot2color'}) && trackingOptions.NCh==2
      trackingOptions.spotsPerNucleus = [1,1];      
      
  elseif ismember(liveExperiment.experimentType,{'inputoutput'}) && trackingOptions.NCh <= 3
    
      % Figure out which channel is the cluster channel
      spotsChannels = liveExperiment.spotChannels;
      inputChannels = liveExperiment.inputChannels;  
      ClusterFilter = ismember(spotsChannels,inputChannels);

      % Cluster channel not limited in spotsPerNucleus        
      trackingOptions.spotsPerNucleus(~ClusterFilter) = 1;

      % Place reasonable limit on link cost for clusters
      trackingOptions.matchCostMax(ClusterFilter) = trackingOptions.matchCostDefault; % this seems to work reasonably well 

      % Place ceiling on number of frames a cluster an go undetected
      % before being "capped"
      trackingOptions.maxUnobservedFrames(ClusterFilter) = ceil(3*60 / Tres); 
      
      % adjust match fraction for clusters
      trackingOptions.targetMatchFrac(ClusterFilter) = 0.9;
  else
      error(['''',liveExperiment.experimentType,''' liveExperiment.experimentType not supported by track04StitchTracks'])
  end