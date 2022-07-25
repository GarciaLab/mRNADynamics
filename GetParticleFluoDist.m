function ParticleFluos = GetParticleFluoDist(CompiledParticles)
%% Calculates all particle schnitz distances for schnitz/particle pairs that are well within the embryo boundaries
if iscell(CompiledParticles)
    CompiledParticles = CompiledParticles{1};
end
ParticleFluos = [];
for p = 1:length(CompiledParticles)
    CurrentParticle = CompiledParticles(p);
    
    ParticleFluos = [ParticleFluos CurrentParticle.Fluo];
    
end