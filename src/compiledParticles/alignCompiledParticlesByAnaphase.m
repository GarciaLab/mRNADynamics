function alignCompiledParticlesByAnaphase(resultsFolder)

load([resultsFolder, filesep, 'CompiledParticles.mat']);
load([resultsFolder, filesep, 'APDivision.mat']);

for ch = 1:length(CompiledParticles)
    for p = 1:length(CompiledParticles{ch})
        [~, apbin] = min(abs(APbinID -  CompiledParticles{ch}(p).MedianAP));
        [~, dvbin] = min(abs(DVbinID -  CompiledParticles{ch}(p).MedianDV));
        divFrames = APDivision(:, apbin);
        actualFrames = CompiledParticles{ch}(p).Frame;
        inds = find(actualFrames(1) > divFrames);
        nc = inds(end);
        CompiledParticles{ch}(p).FramesWRTAnaphase = actualFrames - divFrames(nc);
        CompiledParticles{ch}(p).cycle = nc;
        CompiledParticles{ch}(p).apbin = apbin;
        CompiledParticles{ch}(p).dvbin = dvbin;
    end
end

save([resultsFolder, filesep, 'CompiledParticles.mat'],'CompiledParticles','-append');
