function alignCompiledParticlesByAnaphase(Prefix)



[~,~,DropboxFolder,~, PreProcPath,...
    ~, ~, ~, ~, ~,~] = readMovieDatabase(Prefix);
resultsFolder = [DropboxFolder, filesep, Prefix];

load([resultsFolder, filesep, 'CompiledParticles.mat']);
load([resultsFolder, filesep, 'APDivision.mat']);

for ch = 1:length(CompiledParticles)
    for p = 1:length(CompiledParticles{ch})
        apdif = CompiledParticles{ch}(p).MedianAP - APbinID;
        [~, apbin] = min(apdif(apdif > 0));
        dvdif = CompiledParticles{ch}(p).MedianDV - DVbinID;
        [~, dvbin] = min(dvdif(dvdif > 0));
        divFrames = APDivision(:, apbin);
        if sum(divFrames) == 0
            error('rerun checkdivisiontimes');
        end
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
