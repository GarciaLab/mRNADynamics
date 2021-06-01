function twinParticle = findTwinParticle(CurrentParticle, Particles, CurrentChannelIndex)
if isfield(Particles{CurrentChannelIndex}, 'Schnitz')
    try
    CurrentSchnitz = Particles{CurrentChannelIndex}(CurrentParticle).Schnitz;
    AllParticles = 1:length(Particles{CurrentChannelIndex});
    
    twinParticle = find(([Particles{CurrentChannelIndex}(:).Schnitz] == CurrentSchnitz) & ...
        AllParticles ~= CurrentParticle, 1);
    if isempty(twinParticle)
        twinParticle = [];
    end
    catch
      twinParticle = [];
    end
else
    twinParticle = [];
end
end