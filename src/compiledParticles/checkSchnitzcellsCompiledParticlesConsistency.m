function checkSchnitzcellsCompiledParticlesConsistency(schnitzcells, CompiledParticles)

if iscell(CompiledParticles)
    nSpotChannels = length(CompiledParticles); 
else
    CompiledParticles = {CompiledParticles};
    nSpotChannels = 1;
end

% check for inconsistencies between compiledparticles and schnitzcells
% check that all the schnitz IDs in CompiledParticles are smaller than the
% total number of schnitz in this experiment
for ch = 1:nSpotChannels
    if ~isempty(CompiledParticles{ch}) && isfield(CompiledParticles{ch},'schnitz')
        assert(all([CompiledParticles{ch}.schnitz] <= length(schnitzcells)),'unaccounted schnitz')
    end
end

%if this assert fails, you should run these lines (roughly)
%{
TrackNuclei(prefix,'retrack', 'nWorkers', 1);
trackmRNADynamics(prefix)
CompileParticles(prefix,  'minBinSize', 0, 'MinParticles', 0,'yToManualAlignmentPrompt');
%}

end