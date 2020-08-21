function cptState = updateAuxiliaryStatus(cptState,Status)
  % update auxiliary particle structures
  rawParticleIDs = cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle).idVec;
  rawParticleIDs = unique(rawParticleIDs(~isnan(rawParticleIDs))); % list of indices

  % update Approve flags in auxiliary structures (should probably do
  % this inside CheckParticleTracking eventually
  for r = 1:length(rawParticleIDs)
    cptState.RawParticles{cptState.CurrentChannelIndex}(rawParticleIDs(r)).Approved = Status;
    cptState.HMMParticles{cptState.CurrentChannelIndex}(rawParticleIDs(r)).Approved = Status;
    cptState.SimParticles{cptState.CurrentChannelIndex}(rawParticleIDs(r)).Approved = Status;
  end

  % flag stitch events that correspond to current particle
  linksToUpdate = cellfun(@(x) all(ismember(x,rawParticleIDs)),cptState.ParticleStitchInfo{cptState.CurrentChannelIndex}.linkAdditionIDCell);
  cptState.ParticleStitchInfo{cptState.CurrentChannelIndex}.linkApprovedVec(linksToUpdate) = Status;