function checkSchnitzcellsCompiledParticlesConsistency(schnitzcells, CompiledParticles)

if iscell(CompiledParticles)
    nSpotChannels = length(CompiledParticles); 
else
    CompiledParticles = {CompiledParticles};
    nSpotChannels = 1;
end

%check for inconsistencies between compiledparticles and schnitzcells
for ch = 1:nSpotChannels
    assert( all([CompiledParticles{ch}.schnitz] <= length(schnitzcells)) )
end

%if this assert fails, you should run these lines (roughly)
%{
TrackNuclei(prefix,'retrack', 'nWorkers', 1);
trackmRNADynamics(prefix)
CompileParticles(prefix,  'minBinSize', 0, 'MinParticles', 0,'yToManualAlignmentPrompt');
%}

end