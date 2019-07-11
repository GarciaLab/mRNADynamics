function plotByDorsalConc(Prefix)


[~,~,DropboxFolder,~, PreProcPath,...
    ~, ~, ~, ~, ~,~] = readMovieDatabase(Prefix);
resultsFolder = [DropboxFolder, filesep, Prefix];

load([resultsFolder, filesep, 'CompiledNuclei.mat']);
load([resultsFolder, filesep, 'CompiledParticles.mat']);
load([resultsFolder, filesep, 'APDivision.mat']);

ch = 1;

allRNAs = {[], [], []};
allDorsals = {[], [], []};
ncs = [nc12, nc13, nc14, length(ElapsedTime)];

for p = 1:length(CompiledParticles{ch})
    for nc = 12:14
        if CompiledParticles{ch}(p).cycle == nc
           allRNAs{nc-11} = [allRNAs{nc-11},sum(CompiledParticles{ch}(p).Fluo)];
            n =CompiledParticles{ch}(p).Nucleus
            nucleus = CompiledNuclei(n);
            maxDorsal = max(nucleus.FluoTimeTrace);
            allDorsals{nc-11} = [allDorsals{nc-11}, max(nucleus.FluoTimeTrace)];
        end
    end
end

save([resultsFolder, filesep, 'CompiledParticles.mat'], 'allRNAs', 'allDorsals', '-append');


end