function cptState = updateAuxiliaryStatus(cptState,oldStatus,newStatus)

  % update auxiliary particle structures
  rawParticleIDs = cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).idVec;
  rawParticleIDs = unique(rawParticleIDs(~isnan(rawParticleIDs))); % list of indices
  FragmentIDIndex = [cptState.SimParticles{cptState.CurrentChannelIndex}.FragmentID];
  
  % update Approve flags in auxiliary structures (should probably do
  % this inside CheckParticleTracking eventually
  for r = 1:length(rawParticleIDs)
    cptState.RawParticles{cptState.CurrentChannelIndex}(rawParticleIDs(r)==FragmentIDIndex).Approved = newStatus;
    cptState.HMMParticles{cptState.CurrentChannelIndex}(rawParticleIDs(r)==FragmentIDIndex).Approved = newStatus;
    cptState.SimParticles{cptState.CurrentChannelIndex}(rawParticleIDs(r)==FragmentIDIndex).Approved = newStatus;
  end

  % flag stitch events that correspond to current particle
  linksToUpdate = cellfun(@(x) all(ismember(x,rawParticleIDs)),cptState.ParticleStitchInfo{cptState.CurrentChannelIndex}.linkAdditionIDCell);
  cptState.ParticleStitchInfo{cptState.CurrentChannelIndex}.linkApprovedVec(linksToUpdate) = newStatus;
  if oldStatus == 0
    cptState.ParticleStitchInfo{cptState.CurrentChannelIndex}.reservedFragmentIDs = [cptState.ParticleStitchInfo{cptState.CurrentChannelIndex}.reservedFragmentIDs rawParticleIDs];
  elseif newStatus == 0 
    currIDs = cptState.ParticleStitchInfo{cptState.CurrentChannelIndex}.reservedFragmentIDs;
    cptState.ParticleStitchInfo{cptState.CurrentChannelIndex}.reservedFragmentIDs = currIDs(~ismember(currIDs,rawParticleIDs));
  end