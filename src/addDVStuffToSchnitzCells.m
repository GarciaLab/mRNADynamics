function addDVStuffToSchnitzCells(DataType, varargin)
%%
displayFigures = false;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'displayFigures')
        displayFigures = true;
    end
end
[allData, Prefixes, resultsFolder] = LoadMS2Sets(DataType, 'noCompiledNuclei');

load([resultsFolder,filesep,Prefixes{1},filesep,'FrameInfo.mat'], 'FrameInfo')
imSize = [FrameInfo(1).LinesPerFrame, FrameInfo(1).PixelsPerLine];

if displayFigures
    tileFig = figure();
    tiledlayout('flow');
    holdFig = figure();
end



ch = 1;

fluoFeatures = [];
for e = 1:length(allData)
    
    schnitzcells = allData(e).Particles.schnitzcells;
    
    %to speed this up, remove unnecessary schnitzcells fields
    schnitzcells = removeSchnitzcellsFields(schnitzcells);
    
    CompiledParticles = allData(e).Particles.CompiledParticles;
    DVbinID = allData(1).Particles.DVbinID;
    Ellipses = allData(e).Particles.Ellipses;
    
    for p = 1:length(CompiledParticles{ch})
        schnitzInd = CompiledParticles{ch}(p).schnitz;
        schnitzcells(schnitzInd).compiledParticle = uint16(p);
        if isfield(CompiledParticles{ch}(p), 'dvbin')
            schnitzcells(schnitzInd).dvbin = uint8(CompiledParticles{ch}(p).dvbin);
        end
    end
    
    ncs = [zeros(1,8),allData(e).Particles.nc9, allData(e).Particles.nc10, allData(e).Particles.nc11,...
        allData(e).Particles.nc12, allData(e).Particles.nc13, allData(e).Particles.nc14];
    
    nFrames = length(allData(e).Particles.ElapsedTime);
    
    for s = 1:length(schnitzcells)
        midFrame = ceil(length(schnitzcells(s).frames)/2);
        dif = double(schnitzcells(s).frames(midFrame)) - ncs;
        cycle = max(find(dif>0));
        schnitzcells(s).cycle = uint8(cycle);
    end
    
    
    schnitzcells = filterSchnitz(schnitzcells, imSize);
    
    
    for s = 1:length(schnitzcells)
        schnitzcells(s).FluoTimeTrace = single(ExtractDlFluo(schnitzcells(s).Fluo, .5));
        if schnitzcells(s).cycle ~= 14
            midCycle = floor((ncs(schnitzcells(s).cycle) + ncs(schnitzcells(s).cycle+1))/2);
        else
            midCycle = 1;
        end
        midCycleFrame = find(schnitzcells(s).frames==midCycle);
        schnitzcells(s).FluoFeature = single(schnitzcells(s).FluoTimeTrace(midCycleFrame));
        
        
        
        
        if isempty(schnitzcells(s).FluoFeature)
            schnitzcells(s).FluoFeature = NaN;
        end
        
        if schnitzcells(s).Approved
            
            fluoFeatures = [fluoFeatures, schnitzcells(s).FluoFeature];
            
            if displayFigures && schnitzcells(s).cycle == 12 && ~isnan(schnitzcells(s).FluoFeature)
                
                figure(tileFig)
                nexttile;
                plot(schnitzcells(s).frames, schnitzcells(s).FluoTimeTrace, '-k');
                hold on
                plot(midCycle, schnitzcells(s).FluoFeature, 'ob');
                hold off
                xticks([]);
                yticks([min(schnitzcells(s).FluoTimeTrace), max(schnitzcells(s).FluoTimeTrace)]);
                
                figure(holdFig)
                plot(1:length(schnitzcells(s).FluoTimeTrace), schnitzcells(s).FluoTimeTrace, '-k');
                hold on
                plot(midCycleFrame, schnitzcells(s).FluoFeature, 'ob');
                
                
            end
            
        end
        
    end
    
    if displayFigures
        mkdir([resultsFolder, filesep, DataType]);
        
        saveas(tileFig, [resultsFolder,filesep,DataType, filesep, 'allDorsalTracesTile.png']);
        saveas(tileFig, [resultsFolder,filesep,DataType, filesep, 'allDorsalTracesTile.fig']);
        saveas(tileFig, [resultsFolder,filesep,DataType, filesep, 'allDorsalTracesTile.eps']);
        
        saveas(holdFig, [resultsFolder,filesep,DataType, filesep, 'allDorsalTraces.png']);
        saveas(holdFig, [resultsFolder,filesep,DataType, filesep, 'allDorsalTraces.fig']);
        saveas(holdFig, [resultsFolder,filesep,DataType, filesep, 'allDorsalTraces.eps']);
    end
    
    save([resultsFolder,filesep,Prefixes{e},filesep,Prefixes{e},'_lin.mat'], 'schnitzcells');
    save([resultsFolder,filesep,Prefixes{e},filesep,'Ellipses.mat'], 'Ellipses');
    
end

mkdir([resultsFolder,filesep,DataType]);
% nbins = floor(length(DVbinID) / 2
% nbins = 10;
% dlfluobinwidth = (max(fluoFeatures) - min(fluoFeatures)) / (nbins-1);
% dlfluobins = min(fluoFeatures):dlfluobinwidth:max(fluoFeatures);
dlfluobins = 0:250:3500; %this is appropriate for taking instantaneous dorsal at ~50% through nc12 on the sp8
save([resultsFolder,filesep,DataType,filesep,'dlfluobins.mat'], 'dlfluobins');


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
            schnitzcells(s).dlfluobin = single(dlfluobin);
        else
            schnitzcells(s).dlfluobin = NaN;
        end
    end
    
    
    save([resultsFolder,filesep,Prefixes{e},filesep,Prefixes{e},'_lin.mat'], 'schnitzcells')
    
    for p = 1:length(CompiledParticles{ch})
        schnitzInd = CompiledParticles{ch}(p).schnitz;
        CompiledParticles{ch}(p).dlfluobin = single(schnitzcells(schnitzInd).dlfluobin);
    end
    
    save([resultsFolder,filesep,Prefixes{e},filesep,'CompiledParticles.mat'], 'CompiledParticles', '-append');
    
end

end
