function [dlfluobins, dlfluobincounts] = addDVStuffToSchnitzCells(DataType, varargin)
%%
displayFigures = false;
saveFigures = false;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'displayFigures')
        displayFigures = true;
    elseif strcmpi(varargin{i}, 'saveFigures')
        saveFigures = true;
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
    
    load([resultsFolder,filesep,Prefixes{e},filesep,'FrameInfo.mat'], 'FrameInfo')
    
    ncFrames = [zeros(1,8), allData(e).Particles.nc9, allData(e).Particles.nc10, allData(e).Particles.nc11, allData(e).Particles.nc12, allData(e).Particles.nc13, allData(e).Particles.nc14]; 

    schnitzcells = allData(e).Particles.schnitzcells;
    time = [FrameInfo.Time]/60; %frame times in minutes 
    

    %to speed this up, remove unnecessary schnitzcells fields
    schnitzcells = removeSchnitzcellsFields(schnitzcells);
    
    CompiledParticles = allData(e).Particles.CompiledParticles;
    DVbinID = allData(1).Particles.DVbinID;
    Ellipses = allData(e).Particles.Ellipses;
    
    
    %clear out any existing compiledparticle entries
    for s = 1:length(schnitzcells)
        schnitzcells(s).compiledParticle = [];
    end
        
    
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
    
        schnitzcells = addRelativeTimeToSchnitzcells(schnitzcells, FrameInfo, ncFrames);

    schnitzcells = filterSchnitz(schnitzcells, imSize);
    
    
    for s = 1:length(schnitzcells)
        schnitzcells(s).FluoTimeTrace = single(ExtractDlFluo(schnitzcells(s).Fluo, .5));
        if schnitzcells(s).cycle ~= 14
            midCycle = floor((ncs(schnitzcells(s).cycle) + ncs(schnitzcells(s).cycle+1))/2);
        else
            midCycle = 1;
        end
        midCycleFrame = find(schnitzcells(s).frames==midCycle);
        schnitzcells(s).midCycleFrame = find(schnitzcells(s).frames==midCycle);
        
%         try
            schnitzcells(s).FluoTimeTraceSmooth = smooth(single(ExtractDlFluo(schnitzcells(s).Fluo, .5)), 5);
%         catch
%             schnitzcells(s).FluoTimeTraceSmooth = smooth(single(ExtractDlFluo(schnitzcells(s).Fluo, .5)));
%         end
        if schnitzcells(s).cycle ~= 14
            midCycleSmooth = floor((ncs(schnitzcells(s).cycle) + ncs(schnitzcells(s).cycle+1))/2);
        else
            midCycleSmooth = 1;
        end
        midCycleFrameSmooth = find(schnitzcells(s).frames==midCycleSmooth);
           
        
        %different characterizations of intensity
        schnitzcells(s).fluoMid= single(schnitzcells(s).FluoTimeTrace(midCycleFrame)); %middle fluorescence
        schnitzcells(s).fluoMidSmooth = single(schnitzcells(s).FluoTimeTraceSmooth(midCycleFrameSmooth)); %middle fluorescence smoothed
        
        schnitzcells(s).FluoFeature =  schnitzcells(s).fluoMid; %the preferred intensity
        
        
        
        
        
        
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
                yticks([round(min(schnitzcells(s).FluoTimeTrace)), round(max(schnitzcells(s).FluoTimeTrace))]);
%                 ylim([0, 3000]);
                
                figure(holdFig)
                plot(1:length(schnitzcells(s).FluoTimeTrace), schnitzcells(s).FluoTimeTrace, '-k');
                hold on
                plot(midCycleFram, schnitzcells(s).FluoFeature, 'ob');
%                 ylim([0, 3000]);
                
                
            end
            
        end
        
    end
    
    if displayFigures && saveFigures
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

%%

binDorsal(DataType, false);

end
