function compiledProject = makeCompiledProject(Prefix)

[~, resultsFolder] = getDorsalFolders;

load([resultsFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo')
load([resultsFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'], 'schnitzcells')
load([resultsFolder,filesep,Prefix,filesep,'CompiledParticles.mat'], 'CompiledParticles')

if iscell(CompiledParticles)
    CompiledParticles = CompiledParticles{1};
end



if ~isfield(schnitzcells, 'TimeSinceAnaphase')
    
    ncFrames = getNCFrames(Prefix, resultsFolder);
    schnitzcells = addRelativeTimeToSchnitzcells(...
        schnitzcells, FrameInfo, ncFrames);
    
    CompiledParticles = addRelativeTimeToCompiledParticles(...
        schnitzcells, CompiledParticles, ncFrames, FrameInfo);
end




compiledProject = [];

approvedSchnitzes = find([schnitzcells.Approved]);

n = 0;
for s = approvedSchnitzes
    
        n = n + 1;
        %nuclear compilation
        compiledProject(n).nuclearFrames = schnitzcells(s).frames;
        compiledProject(n).cycle = schnitzcells(s).cycle;
        compiledProject(n).compiledParticle = schnitzcells(s).compiledParticle;
        compiledProject(n).dorsalFluoBin = schnitzcells(s).dorsalFluoBin;
        compiledProject(n).dorsalFluoFeature = schnitzcells(s).FluoFeature;
        compiledProject(n).dorsalFluoTimeTrace= schnitzcells(s).FluoTimeTrace;
        compiledProject(n).nuclearTimeSinceAnaphase = schnitzcells(s).timeSinceAnaphase;
        
        
        %particle compilation
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


save([resultsFolder,filesep,Prefix,filesep,'compiledProject.mat'], 'compiledProject');


end