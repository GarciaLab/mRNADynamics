function addDVStuffToSchnitzCells(DataType)
%%
[allData, Prefixes, resultsFolder] = LoadMS2Sets(DataType);

ch = 1;

fluoFeatures = [];
for e = 1:length(allData)
    
    schnitzcells = allData(e).Particles.schnitzcells;
    CompiledParticles = allData(e).Particles.CompiledParticles;
    DVbinID = allData(1).Particles.DVbinID;
    
    for p = 1:length(CompiledParticles)
        schnitzInd = CompiledParticles{ch}(p).schnitz;
        schnitzcells(schnitzInd).compiledParticle = p;
        schnitzcells(schnitzInd).cycle = CompiledParticles{ch}(p).cycle;
        schnitzcells(schnitzInd).dvbin = CompiledParticles{ch}(p).dvbin;
    end
    
    ncs = [zeros(1,8),allData(e).Particles.nc9, allData(e).Particles.nc10, allData(e).Particles.nc11,...
        allData(e).Particles.nc12, allData(e).Particles.nc13, allData(e).Particles.nc14];
    
    for s = 1:length(schnitzcells)
        dif = schnitzcells(s).frames(1) - ncs;
        cycle = max(find(dif>0));
        schnitzcells(s).cycle = cycle;           
    end
    
    for s = 1:length(schnitzcells)
        schnitzcells(s).FluoTimeTrace = ExtractDlFluo(schnitzcells(s).Fluo, .5);
        if schnitzcells(s).cycle ~= 14
            midCycle = floor((ncs(schnitzcells(s).cycle) + ncs(schnitzcells(s).cycle+1))/2);
        else
            midCycle = 1;
        end
        midCycleFrame = find(schnitzcells(s).frames==midCycle);
        schnitzcells(s).FluoFeature = schnitzcells(s).FluoTimeTrace(midCycleFrame);
        if isempty(schnitzcells(s).FluoFeature)
            schnitzcells(s).FluoFeature = NaN;
        end
        fluoFeatures = [fluoFeatures, schnitzcells(s).FluoFeature];
    end
  
    save([resultsFolder,filesep,Prefixes{e},filesep,Prefixes{e},'_lin.mat'], 'schnitzcells')
    
    
end

nbins = length(DVbinID);
dlfluobinwidth = (max(fluoFeatures) - min(fluoFeatures)) / (nbins-1);
dlfluobins = min(fluoFeatures):dlfluobinwidth:max(fluoFeatures);

%%

[allData, Prefixes, resultsFolder] = LoadMS2Sets(DataType);
load([resultsFolder,filesep,DataType,filesep,'dlfluobins.mat'], 'dlfluobins');

for e = 1:length(allData)
    
    schnitzcells = allData(e).Particles.schnitzcells;
    CompiledParticles = allData(e).Particles.CompiledParticles;
    
    for s = 1:length(schnitzcells)
        dif = schnitzcells(s).FluoFeature - dlfluobins;
        [~,dlfluobin] = min(dif(dif>0));
        if ~isempty(dlfluobin)
            schnitzcells(s).dlfluobin = dlfluobin;
        else
            schnitzcells(s).dlfluobin = NaN;
        end
    end
    

    save([resultsFolder,filesep,Prefixes{e},filesep,Prefixes{e},'_lin.mat'], 'schnitzcells')
    
    for p = 1:length(CompiledParticles{ch})
        schnitzInd = CompiledParticles{ch}(p).schnitz;
        CompiledParticles{ch}(p).dlfluobin = schnitzcells(schnitzInd).dlfluobin;
    end
    
    save([resultsFolder,filesep,Prefixes{e},filesep,'CompiledParticles.mat'], 'CompiledParticles', '-append');
    
end

end
