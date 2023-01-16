function MaxTraceLengths = GetParticleTraceLengthDist(CompiledParticles)
%% Calculates all particle schnitz distances for schnitz/particle pairs that are well within the embryo boundaries
if iscell(CompiledParticles)
    CompiledParticles = CompiledParticles{1};
end
TraceLengths = [];
cycles = [];
for p = 1:length(CompiledParticles)
    CurrentParticle = CompiledParticles(p);
    if ~isempty(CurrentParticle.cycle) & CurrentParticle.cycle > 0 & isfield(CurrentParticle.schnitzcell, 'FlaggingInfo')
        TraceLengths(end+1) = length(CurrentParticle.schnitzcell.FlaggingInfo.AllSchnitzFrames(CurrentParticle.schnitzcell.FlaggingInfo.SchnitzPresent));
        cycles(end+1) = CurrentParticle.cycle;
    end
    
end
MaxTraceLengths = NaN(1, 6);
for NC = 9:14
    MaxTraceLengths(NC-8) = prctile(TraceLengths(cycles == NC), 95);
end