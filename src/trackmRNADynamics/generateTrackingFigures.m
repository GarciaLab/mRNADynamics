
function [ParticlesFig, particlesAxes, NucleiFig, nucAxes] = generateTrackingFigures(app, UseHistone)

ParticlesFig=[]; particlesAxes=[]; NucleiFig=[];nucAxes = [];

if isempty(app)
    
    ParticlesFig = figure;
    particlesAxes = axes(ParticlesFig);
    
    if UseHistone
        NucleiFig = figure;
        set(NucleiFig, 'units', 'normalized', 'position', [0.65, .5, .2, .2])
        nucAxes = axes(NucleiFig);
    end
    
end

end