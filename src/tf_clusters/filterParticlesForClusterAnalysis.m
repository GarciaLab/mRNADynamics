function filteredParticles = filterParticlesForClusterAnalysis(...
                        Particles, ncStart, minTraceLen)

filteredParticles = struct();
counter = 1;
frame = [Particles.Frame];
traceLen = numel(frame);
traces = frame(frame>=ncStart && frame)
for i = 1:length(Particles)
    traceStart = Particles(i).Frame(1);
    traceLen = length(Particles(i).Frame);
    
    if (traceStart <= ncStart) || (traceLen < minTraceLen)
        filteredParticles(i) = true;
    end
end

if isempty(filteredParticles)
    error('No particles met your filtering criteria. Please adjust inputs and try again.')
end

