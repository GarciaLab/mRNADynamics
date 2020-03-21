function compiledProject = makeCompiledProject(Prefix)

[~, resultsFolder] = getDorsalFolders;

ncFrames = getNCFrames(Prefix, resultsFolder);

load([resultsFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo')

ncFrames(ncFrames==0) = 1;
time = [FrameInfo.Time]/60; %frame times in minutes
ncFrames(isnan(ncFrames)) = [];
ncTimes = time(ncFrames);
if length(ncTimes) < 14
    ncTimes(end+1:14) = nan;
end

load([resultsFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'], 'schnitzcells')
load([resultsFolder,filesep,Prefix,filesep,'CompiledParticles.mat'], 'CompiledParticles')

if iscell(CompiledParticles)
    CompiledParticles = CompiledParticles{1};
end

for s = 1:length(schnitzcells)
    schnitzcells(s).timeSinceAnaphase = time(schnitzcells(s).frames) - ncTimes(schnitzcells(s).cycle);
end

clear s;

for p = 1:length(CompiledParticles)
    s = CompiledParticles(p).Nucleus;
    CompiledParticles(p).timeSinceAnaphase = time(CompiledParticles(p).Frame) - ncTimes(schnitzcells(s).cycle);
    if CompiledParticles(p).timeSinceAnaphase == 0 & schnitzcells(s).cycle == 12
        'stop'
    end
end


compiledProject = [];
n = 0;
for s = 1:length(schnitzcells)
    if schnitzcells(s).Approved
        n = n + 1;
        %nuclear stuff
        compiledProject(n).nuclearFrames = schnitzcells(s).frames;
        compiledProject(n).cycle = schnitzcells(s).cycle;
        compiledProject(n).compiledParticle = schnitzcells(s).compiledParticle;
        compiledProject(n).dorsalFluoBin = schnitzcells(s).dorsalFluoBin;
        compiledProject(n).dorsalFluoFeature = schnitzcells(s).FluoFeature;
        compiledProject(n).dorsalFluoTimeTrace= schnitzcells(s).FluoTimeTrace;
        compiledProject(n).nuclearTimeSinceAnaphase = schnitzcells(s).timeSinceAnaphase;
        
        
        %particle stuff
        p = schnitzcells(s).compiledParticle;
        if ~isempty(p)
            if CompiledParticles(p).Approved
                compiledProject(n).particleFrames = CompiledParticles(p).Frame;
                compiledProject(n).particleTimeSinceAnaphase = CompiledParticles(p).timeSinceAnaphase;
                
                
                
                compiledProject(n).particleFluo= CompiledParticles(p).Fluo;
                compiledProject(n).particleFluo3Slice = CompiledParticles(p).Fluo3;
                compiledProject(n).particleOffset = CompiledParticles(p).Off;
                compiledProject(n).particleFluo3D = CompiledParticles(p).Fluo3DRaw;
                
                fluo = compiledProject(n).particleFluo3D;
                tau = compiledProject(n).particleTimeSinceAnaphase;
                
                %some preliminary additional statistics about the spot
                compiledProject(n).particleDuration = max(tau) - min(tau);
                compiledProject(n).particleFluo95= fluo(fluo>=prctile(fluo,95));
                compiledProject(n).particleTimeOn = min(tau);
                
                if length(tau) > 1
                    compiledProject(n).particleAccumulatedFluo = trapz(tau, fluo, 2);
                else
                    compiledProject(n).particleAccumulatedFluo = fluo;
                end
                
            else
                compiledProject(n).particleFrames = [];
                compiledProject(n).particleFluo= [];
                compiledProject(n).particleFluo3Slice = [];
                compiledProject(n).particleOffset = [];
                compiledProject(n).particleFluo3D = [];
                compiledProject(n).particleDuration = [];
                compiledProject(n).particleFluo95= [];
                compiledProject(n).particleTimeOn = [];
                compiledProject(n).particleAccumulatedFluo = [];
            end
        end
    end
end


save([resultsFolder,filesep,Prefix,filesep,'compiledProject.mat'], 'compiledProject');


end