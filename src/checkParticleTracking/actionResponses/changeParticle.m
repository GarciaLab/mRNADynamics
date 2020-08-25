function cptState =  changeParticle(ParticleNum, cptState)
  
    Particles = cptState.Particles;
    numParticles = cptState.numParticles;
    CurrentChannel = cptState.CurrentChannelIndex;

    cptState.CurrentParticle = min(max(ParticleNum, 1), numParticles);
    cptState.ManualZFlag = 0;
    
    % check for QC flags
    allFlagVec = false(1,length(Particles{CurrentChannel}(cptState.CurrentParticle).Frame));    
    allUrgentFlagVec = false(1,length(Particles{CurrentChannel}(cptState.CurrentParticle).Frame));    
    for v = 1:length(cptState.qcFlagFields)
      allUrgentFlagVec = allUrgentFlagVec | Particles{CurrentChannel}(cptState.CurrentParticle).(cptState.qcFlagFields{v})>1;
      allFlagVec = allFlagVec | Particles{CurrentChannel}(cptState.CurrentParticle).(cptState.qcFlagFields{v})>0;      
    end
    if any(allUrgentFlagVec)
      cptState.CurrentFrame = Particles{CurrentChannel}(cptState.CurrentParticle).Frame(find(allUrgentFlagVec,1));
    elseif any(allFlagVec)
      cptState.CurrentFrame = Particles{CurrentChannel}(cptState.CurrentParticle).Frame(find(allFlagVec,1));
    else
      cptState.CurrentFrame = Particles{CurrentChannel}(cptState.CurrentParticle).Frame(1);
    end
end

