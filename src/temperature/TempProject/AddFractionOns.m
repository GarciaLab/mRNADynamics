function this = AddFractionOns(this)
% author: G. Martini
% date created: 1/19/21
% date last modified: 1/19/21
NumExperiments = length(this.ExperimentPrefixes);
APResolution = this.Experiments{1}.APResolution;
NumAPbins = uint16(1/APResolution)+1;
APboundaries = 0:APResolution:1;
APlowers = APboundaries(1:end-1);
APuppers = APboundaries(2:end);


this.SchnitzCounts = zeros(NumExperiments, NumAPbins, 6);
this.FractionOns = NaN(NumExperiments, NumAPbins, 6);
for SetIndex = 1:length(this.ExperimentPrefixes)
    %disp(num2str(SetIndex))
    if ~ismember(SetIndex, this.ProcessedExperiments)
        continue
    end
    schnitzcells = getSchnitzcells(this.Experiments{SetIndex});
    Particles = getParticles(this.Experiments{SetIndex});
    if iscell(Particles) & (length(Particles) == 1)
        Particles = Particles{1};
    end
    load([this.Experiments{SetIndex}.resultsFolder, 'CompiledParticles.mat'], 'CompiledParticles')
    if iscell(CompiledParticles) & (length(CompiledParticles) == 1)
        CompiledParticles = CompiledParticles{1};
    end
    
    CompiledParticleSchnitzMap = [CompiledParticles(:).Nucleus];
    if (length(unique(CompiledParticleSchnitzMap)) ~= length(CompiledParticleSchnitzMap)) | ...
            (length(CompiledParticleSchnitzMap) ~= length(CompiledParticles)) 
        %disp(num2str(SetIndex))
        %warning('Compiled Particle mapping to schnitz cells either contains duplicates or unmapped particles.')
        CompiledParticles = ...
            correctParticleMappingsToSchnitzCells(this.ExperimentPrefixes{SetIndex});
        if iscell(CompiledParticles) & (length(CompiledParticles) == 1)
            CompiledParticles = CompiledParticles{1};
        end
        CompiledParticleSchnitzMap = [CompiledParticles(:).Nucleus];
    
    end
    FrameInfo = getFrameInfo(this.Experiments{SetIndex});
    PixelSize = this.Experiments{SetIndex}.pixelSize_um;
    nucleusDiameters = zeros(1, 6);
    for nc=9:14
        nucleusDiameters(nc-8) = getDefaultParameters(FrameInfo,[['d', num2str(nc)]])/PixelSize; % in pixels
    end
    xDim = this.Experiments{SetIndex}.xDim;
    yDim = this.Experiments{SetIndex}.yDim;
    zDim = this.Experiments{SetIndex}.zDim;
    
    SchnitzApproved = zeros(1, length(schnitzcells));
    HasAllFrames = zeros(1, length(schnitzcells));
    SchnitzInBoundsAllFrames = zeros(1, length(schnitzcells));
    SchnitzCycles = zeros(1, length(schnitzcells));
    SchnitzAPs = zeros(1, length(schnitzcells));
    SchnitzHasValidParticle = zeros(1, length(schnitzcells));
    
    
    
    
    for sc_idx=1:length(schnitzcells)
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
            CurrentCP = CompiledParticles(cp_idx);
            NValidFrames = sum(CurrentCP.ManualFrameApproved(1:length(CurrentCP.Frame)));
            if CurrentCP.ManualApproved & (NValidFrames >= this.MinimumTimePoints)
                SchnitzHasValidParticle(sc_idx) = 1;
            end
        end
    end
    for NC = 9:14
        for APidx = 1:NumAPbins-1
            this.SchnitzCounts(SetIndex,APidx, NC-8)  = sum((SchnitzApproved == 1) &...
                (SchnitzInBoundsAllFrames == 1) & (HasAllFrames == 1) & ...
                (SchnitzCycles == NC) & (SchnitzAPs >= APlowers(APidx)) & (SchnitzAPs < APuppers(APidx)));
            this.FractionOns(SetIndex,APidx, NC-8)  = sum((SchnitzApproved == 1) & ...
                (SchnitzInBoundsAllFrames == 1) & (HasAllFrames == 1) & ...
                (SchnitzCycles == NC) & (SchnitzAPs >= APlowers(APidx)) & (SchnitzAPs < APuppers(APidx)) & ...
                (SchnitzHasValidParticle == 1))/this.SchnitzCounts(SetIndex,APidx, NC-8);
        end
    end
    
end
end