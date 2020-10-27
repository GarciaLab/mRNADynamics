function matchCostVec = determineMatchOptions(Prefix,useHistone,matchCostMax)

  liveExperiment = LiveExperiment(Prefix);
  ExperimentType = liveExperiment.experimentType;
  NCh = length(liveExperiment.spotChannels);
  % simplify the logic a bit
  traceCost = matchCostMax;
  if useHistone
    traceCost = realmax;
  end
  
  if ismember(ExperimentType,{'1spot'})||ismember(ExperimentType,{'2spot'})||ismember(ExperimentType,{'2spot2color'})  
    matchCostVec = repelem(traceCost,NCh);
  elseif ismember(ExperimentType,{'inputoutput'}) 
    % Do we have TF clusters in the input channel?
    % No clusters:
    if NCh == 1
        matchCostVec = traceCost;
    % Clusters:
    elseif NCh == 2
        % Figure out which channel is the cluster channel
        spotsChannels = liveExperiment.spotChannels;
        inputChannels = liveExperiment.inputChannels;
        % Cluster channel should still have a sigma limit on fragment
        % matching
        matchCostVec(spotsChannels == intersect(inputChannels, spotsChannels)) = matchCostMax;
        matchCostvec(spotsChannels ~= intersect(inputChannels, spotsChannels)) = traceCost;
    else
        error('No spot channels, or too many, detected.')
    end
  else
      error(['''',ExperimentType,''' ExperimentType not supported by track04StitchTracks'])
  end