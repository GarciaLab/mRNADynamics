function cptState = removeSpotFromAuxParticles(cptState,Prefix)
  CC = cptState.CurrentChannelIndex;
  CP = cptState.CurrentParticle;
  CF = cptState.CurrentFrame;
  reservedFragmentIDs = cptState.ParticleStitchInfo{CC}.reservedFragmentIDs;    
  
  % remove this spot from all auxiliary particle structure... 
  FragmentIDVec = [cptState.SimParticles{CC}.FragmentID];
  FragmentID = cptState.Particles{CC}(CP).idVec(CF);
  FragmentFilter = FragmentIDVec==FragmentID;    
  overlapFilter = ismember(cptState.SimParticles{CC}(FragmentFilter).Frame,CF);    
  
  % reset Particle's id vector  
  cptState.Particles{CC}(CP).idVec(CF) = NaN;
  
  if all(overlapFilter)
    % Case 1: if the full fragment falls inside region, then remove entirely      
    cptState.SimParticles{CC} = cptState.SimParticles{CC}(updateIndex);       

    % other particle structures
    cptState.RawParticles{CC} = cptState.RawParticles{CC}(~FragmentFilter);
    cptState.HMMParticles{CC} = cptState.RawParticles{CC}(~FragmentFilter);

  elseif sum(~overlapFilter)==1 || max(diff(find(~overlapFilter))) == 1
    % case 2: if a contiguous fragment remains, update entry            
    cptState.RawParticles{CC}(FragmentFilter) = ...
      updateVectorFields(cptState.RawParticles{CC}(FragmentFilter),~overlapFilter);

    % now we need to re-generate motion model and path predictions
    globalMotionModel = getGlobalMotionModel(LiveExperiment(Prefix));
    
    [HMMTemp, ~] = track02TrainGHMM(...
      {cptState.RawParticles{CC}(FragmentFilter)}, globalMotionModel, false);

    SimTemp = track03PredictParticlePaths(HMMTemp, cptState.FrameInfo, false);

    % update structures
    cptState.HMMParticles{CC}(FragmentFilter) = HMMTemp{1};
    cptState.SimParticles{CC}(FragmentFilter) = SimTemp{1};
  else            
    % case 3: in the unlikely event that a hole is created in the 
    % middle of a contiguous fragment, then we must create a new
    % entry for the trailing fragment
    tempIDVec = bwlabel(~overlapFilter);
    tempIDIndex = unique(tempIDVec(tempIDVec~=0));
    newIndices = [find(FragmentFilter) length(FragmentIDVec)+(1:max(tempIDVec-1))];
    RawOrig = cptState.RawParticles{CC}(FragmentFilter);
    for id = 1:length(tempIDIndex)
      subFilter = tempIDVec==tempIDIndex(id);
      % first update raw particles               
      cptState.RawParticles{CC}(newIndices(id)) = ...
                            updateVectorFields(RawOrig,subFilter);

      % change fragment ID if necessary
      ptID = FragmentID;
      if id > 1
        ptID = nanmax(FragmentIDVec)+1;
        while ismember(ptID,reservedFragmentIDs)
          ptID = ptID + 1;
        end        
      end
      % then retrain motion model
      globalMotionModel = getGlobalMotionModel(LiveExperiment(Prefix));
      [HMMTemp, ~] = track02TrainGHMM(...
        {cptState.RawParticles{CC}(newIndices(id))}, globalMotionModel, false);
      % then make path predictions
      SimTemp = track03PredictParticlePaths(HMMTemp, cptState.FrameInfo, false);
      % update structures
      cptState.HMMParticles{CC}(newIndices(id)) = HMMTemp{1};
      cptState.SimParticles{CC}(newIndices(id)) = SimTemp{1};
      % update Particles
      cptState.Particles{CC}(CP).idVec(SimTemp{1}.Frame) = ptID;
    end
  end
  