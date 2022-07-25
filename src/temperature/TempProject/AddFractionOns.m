function ltm = AddFractionOns(ltm)
% author: G. Martini
% date created: 1/19/21
% date last modified: 1/19/21
NumExperiments = length(ltm.ExperimentPrefixes);
APResolution = ltm.Experiments{1}.APResolution;
NumAPbins = uint16(1/APResolution)+1;
APboundaries = 0:APResolution:1;
APlowers = APboundaries(1:end-1);
APuppers = APboundaries(2:end);


ltm.SchnitzCounts = zeros(NumExperiments, NumAPbins, 6);
ltm.FractionOns = NaN(NumExperiments, NumAPbins, 6);
disp('Adding Fraction Ons.')
for SetIndex = 1:length(ltm.ExperimentPrefixes)
    %disp(num2str(SetIndex))
    if ~ismember(SetIndex, ltm.ProcessedExperiments)
        continue
    end
    schnitzcells = getSchnitzcells(ltm.Experiments{SetIndex});
    Particles = getParticles(ltm.Experiments{SetIndex});
    if iscell(Particles) & (length(Particles) == 1)
        Particles = Particles{1};
    end
    load([ltm.Experiments{SetIndex}.resultsFolder, 'CompiledParticles.mat'], 'CompiledParticles')
    if iscell(CompiledParticles) & (length(CompiledParticles) == 1)
        CompiledParticles = CompiledParticles{1};
    end
    
    CompiledParticleSchnitzMap = [CompiledParticles(:).Nucleus];
    if (length(unique(CompiledParticleSchnitzMap)) ~= length(CompiledParticleSchnitzMap)) | ...
            (length(CompiledParticleSchnitzMap) ~= length(CompiledParticles))
        %disp(num2str(SetIndex))
        %warning('Compiled Particle mapping to schnitz cells either contains duplicates or unmapped particles.')
        CompiledParticles = ...
            correctParticleMappingsToSchnitzCells(ltm.ExperimentPrefixes{SetIndex});
        if iscell(CompiledParticles) & (length(CompiledParticles) == 1)
            CompiledParticles = CompiledParticles{1};
        end
        CompiledParticleSchnitzMap = [CompiledParticles(:).Nucleus];
        
    end
    FrameInfo = getFrameInfo(ltm.Experiments{SetIndex});
    PixelSize = ltm.Experiments{SetIndex}.pixelSize_um;
    nucleusDiameters = zeros(1, 6);
    for nc=9:14
        nucleusDiameters(nc-8) = getDefaultParameters(FrameInfo,[['d', num2str(nc)]])/PixelSize; % in pixels
    end
    xDim = ltm.Experiments{SetIndex}.xDim;
    yDim = ltm.Experiments{SetIndex}.yDim;
    zDim = ltm.Experiments{SetIndex}.zDim;
    
    SchnitzApproved = zeros(1, length(schnitzcells));
    HasAllFrames = zeros(1, length(schnitzcells));
    SchnitzInBoundsAllFrames = zeros(1, length(schnitzcells));
    SchnitzCycles = zeros(1, length(schnitzcells));
    SchnitzAPs = zeros(1, length(schnitzcells));
    SchnitzHasValidParticle = zeros(1, length(schnitzcells));
    
    
    
    
    for sc_idx=1:length(schnitzcells)
        if ~isfield(schnitzcells(sc_idx), 'Approved')
            schnitzcells(sc_idx).Approved = 0;
        end
        if (schnitzcells(sc_idx).Approved > 0)  & (schnitzcells(sc_idx).Flag ~= 6)
            SchnitzApproved(sc_idx) = 1;
        end
        if schnitzcells(sc_idx).VelocityInfo.SchnitzHasAllFrames
            HasAllFrames(sc_idx) =  1;
            nc_idx = schnitzcells(sc_idx).cycle-8;
            sc_d = nucleusDiameters(nc_idx);
            IsInBounds = all((schnitzcells(sc_idx).cenx >= sc_d/2) & ...
                (schnitzcells(sc_idx).cenx <= xDim-sc_d/2) & ...
                (schnitzcells(sc_idx).ceny >= sc_d/2) & ...
                (schnitzcells(sc_idx).ceny <= yDim-sc_d/2));
            if IsInBounds
                SchnitzInBoundsAllFrames(sc_idx) = 1;
            end
        end
        
        if ~isempty(schnitzcells(sc_idx).cycle)
            SchnitzCycles(sc_idx) =  schnitzcells(sc_idx).cycle;
        end
        SchnitzAPs(sc_idx) = mean(schnitzcells(sc_idx).APpos);
        cp_idx = find(CompiledParticleSchnitzMap == sc_idx);
        if ~isempty(cp_idx)
            for cp2_idx = 1:length(cp_idx)
                cp2 = cp_idx(cp2_idx);
                CurrentCP = CompiledParticles(cp2);
                NValidFrames = sum(CurrentCP.ManualFrameApproved(1:length(CurrentCP.Frame)));
                if CurrentCP.ManualApproved & (NValidFrames >= ltm.MinimumTimePoints)
                    SchnitzHasValidParticle(sc_idx) = 1;
                end
            end
        end
    end
    for NC = 9:14
        for APidx = 1:NumAPbins-1
            ltm.SchnitzCounts(SetIndex,APidx, NC-8)  = sum((SchnitzApproved == 1) &...
                (SchnitzInBoundsAllFrames == 1) & (HasAllFrames == 1) & ...
                (SchnitzCycles == NC) & (SchnitzAPs >= APlowers(APidx)) & (SchnitzAPs < APuppers(APidx)));
            ltm.FractionOns(SetIndex,APidx, NC-8)  = sum((SchnitzApproved == 1) & ...
                (SchnitzInBoundsAllFrames == 1) & (HasAllFrames == 1) & ...
                (SchnitzCycles == NC) & (SchnitzAPs >= APlowers(APidx)) & (SchnitzAPs < APuppers(APidx)) & ...
                (SchnitzHasValidParticle == 1))/ltm.SchnitzCounts(SetIndex,APidx, NC-8);
        end
    end
    
end

end