function alignCompiledParticlesByAnaphase(Prefix)



[~,~,DropboxFolder,~, PreProcPath,...
    ~, ~, ~, ~, ~,~] = readMovieDatabase(Prefix);
resultsFolder = [DropboxFolder, filesep, Prefix];

fullEmbryoExists =  exist([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'], 'file');

load([resultsFolder, filesep, 'CompiledParticles.mat']);
if fullEmbryoExists
    try
        load([resultsFolder, filesep, 'APDivision.mat']);
    catch
        warning('apdivision not found. making it now.')
        CheckDivisionTimes(Prefix, 'lazy');
        load([resultsFolder, filesep, 'APDivision.mat']);
    end
end


for ch = 1:length(CompiledParticles)
    for p = 1:length(CompiledParticles{ch})
        if fullEmbryoExists
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
        
        try
            CompiledParticles{ch}(p).cycle = nc;
        catch
            CompiledParticles{ch}(p).cycle = CompiledParticles{ch}(p).nc;
        end
        
    end
end



save([resultsFolder, filesep, 'CompiledParticles.mat'],'CompiledParticles','-append');
