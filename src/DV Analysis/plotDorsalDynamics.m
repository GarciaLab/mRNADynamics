function figs = plotDorsalDynamics(DataType, varargin)

figsFlag = false;
col = 'b';

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'figs')
        figs = varargin{i+1};
        figsFlag = true;
        col = 'r';
    end
end

[allData, Prefixes, resultsFolder] = LoadMS2Sets(DataType, 'noCompiledNuclei');

load([resultsFolder,filesep,DataType,filesep,'dlfluobins.mat'], 'dlfluobins');


nc=12;
ch = 1;

for bin = 1:12 %just look at the first 12 bins
    
    if figsFlag
        figure(figs{bin})
    else
        figure();
        figs{bin} = gcf;
    end
    
    for e = 1:length(allData)
        
        load([resultsFolder,filesep,Prefixes{e},filesep,'FrameInfo.mat'], 'FrameInfo')
           
  
        
        CompiledParticles = allData(e).Particles.CompiledParticles;
        schnitzcells = allData(e).Particles.schnitzcells;
        
        ncFrames = [zeros(1,8),allData(e).Particles.nc9, allData(e).Particles.nc10, allData(e).Particles.nc11,...
            allData(e).Particles.nc12, allData(e).Particles.nc13, allData(e).Particles.nc14];
        
        ncFrames(ncFrames==0) = 1;
        time = [FrameInfo.Time]/60; %frame times in minutes 
        ncTimes = time(ncFrames);
        
        if ~isfield(schnitzcells, 'timeSinceAnaphase')
            schnitzcells = addRelativeTimeToSchnitzcells(schnitzcells, FrameInfo, ncFrames);
        end
        
        for s = 1:length(schnitzcells)
            
            if schnitzcells(s).cycle ~= 14
                midCycleFrame = floor((ncFrames(schnitzcells(s).cycle) + ncFrames(schnitzcells(s).cycle+1))/2);
            else
                midCycleFrame = 1;
            end
            
            midCycleFrameIndex = find(schnitzcells(s).frames==midCycleFrame);
            
            if schnitzcells(s).cycle == nc & schnitzcells(s).dlfluobin == bin...
                    & schnitzcells(s).Approved

                yyaxis left
                plot(schnitzcells(s).timeSinceAnaphase, schnitzcells(s).FluoTimeTrace, ['-', col]);
                hold on
                plot(schnitzcells(s).timeSinceAnaphase(midCycleFrameIndex), schnitzcells(s).FluoFeature, ['o', col]);
                
                p = schnitzcells(s).compiledParticle;
                if ~isempty(p) & p < length(CompiledParticles{ch}) & CompiledParticles{ch}(p).Nucleus == s
                    yyaxis right
                    pTime = time(CompiledParticles{ch}(p).Frame) - ncTimes(schnitzcells(s).cycle);
                    pFluo = CompiledParticles{ch}(p).Fluo3DRaw'; 
                
                    plot(pTime, pFluo, ['-.', col]);
                    hold on
                end
                
            end
            
        end
        
    end
    
    ylim([0, 3000]);
    title(num2str(bin));
    
end


end