function [ CompiledParticles, schnitzcells, Ellipses, Particles] = ...
    correctParticleMappingsToSchnitzCells(Prefix,schnitzcells, Ellipses, Particles, CompiledParticles)

%% DESCRIPTION this function corrects mapping errors between schnitzcells, Ellipses, Particles and CompiledParticles
liveExperiment = LiveExperiment(Prefix);
if ~exist('schnitzcells', 'var')
    schnitzcells = getSchnitzcells(liveExperiment);
end
if ~exist('Ellipses', 'var')
    Ellipses = getEllipses(liveExperiment);
end
if ~exist('Particles', 'var')
    Particles = getParticles(liveExperiment);
end
if ~exist('CompiledParticles', 'var')
    load([liveExperiment.resultsFolder, 'CompiledParticles.mat'], 'CompiledParticles')
end

[schnitzcells, Ellipses] = correctSchnitzCellErrors(schnitzcells, Ellipses, Prefix);
for ChN=liveExperiment.spotChannels
    for p = 1:length(Particles{ChN})
        if isempty(Particles{ChN}(p).Nucleus)
            continue
        end
        ClosestSchnitzCellsInFrames = zeros(1, length(Particles{ChN}(p).Frame));
        for frame_index = 1:length(Particles{ChN}(p).Frame)
            spotPos = [Particles{ChN}(p).xPos(frame_index), Particles{ChN}(p).yPos(frame_index)];
            EllipsesPos = [Ellipses{Particles{ChN}(p).Frame(frame_index)}(:,1), Ellipses{Particles{ChN}(p).Frame(frame_index)}(:,2)];
            Distances = vecnorm(spotPos - EllipsesPos, 2, 2).';
            [~, minimumDistanceIndex] = min(Distances);
            if ~isempty(minimumDistanceIndex)
                ClosestSchnitzCellsInFrames(frame_index) = Ellipses{Particles{ChN}(p).Frame(frame_index)}(minimumDistanceIndex,9);
            end
        end
        SchnitzCellIndexForParticle = mode(ClosestSchnitzCellsInFrames);
        if ~isempty(SchnitzCellIndexForParticle)
            Particles{ChN}(p).Nucleus = SchnitzCellIndexForParticle;
        end
    end
    if exist('CompiledParticles', 'var')
        
        
        for i = 1:length(CompiledParticles{ChN})
            
            if isempty(CompiledParticles{ChN}(i).Nucleus)
                if isfield(CompiledParticles, 'schnitzcell') 
                    CompiledParticles{ChN}(i).schnitzcell = {};
                end
                continue
            end
            ClosestSchnitzCellsInFrames = zeros(1, length(CompiledParticles{ChN}(i).Frame));
            for frame_index = 1:length(CompiledParticles{ChN}(i).Frame)
                spotPos = [CompiledParticles{ChN}(i).xPos(frame_index), CompiledParticles{ChN}(i).yPos(frame_index)];
                EllipsesPos = [Ellipses{CompiledParticles{ChN}(i).Frame(frame_index)}(:,1), Ellipses{CompiledParticles{ChN}(i).Frame(frame_index)}(:,2)];
                Distances = vecnorm(spotPos - EllipsesPos, 2, 2).';
                [~, minimumDistanceIndex] = min(Distances);
                if ~isempty(minimumDistanceIndex)
                    ClosestSchnitzCellsInFrames(frame_index) = Ellipses{CompiledParticles{ChN}(i).Frame(frame_index)}(minimumDistanceIndex,9);
                end
            end
            SchnitzCellIndexForParticle = mode(ClosestSchnitzCellsInFrames);
            if ~isempty(SchnitzCellIndexForParticle)
                CompiledParticles{ChN}(i).Nucleus = SchnitzCellIndexForParticle;
                if isfield(CompiledParticles, 'schnitzcell') 
                    CompiledParticles{ChN}(i).schnitzcell = schnitzcells(SchnitzCellIndexForParticle);
                end
            else
                if isfield(CompiledParticles, 'schnitzcell') 
                    CompiledParticles{ChN}(i).schnitzcell = {};
                end
            end
        end
        CompiledParticleSchnitzMap = [CompiledParticles{ChN}(:).Nucleus];
        if length(CompiledParticleSchnitzMap) < length(CompiledParticles{ChN})
            cp_indices = 1:length(CompiledParticles{ChN});
            for i = 1:length(cp_indices)
                cp_idx = cp_indices(i);
                if isempty(CompiledParticles{ChN}(cp_idx).Nucleus)
                    CompiledParticles{ChN}(cp_idx) = [];
                    cp_indices = cp_indice -1;
                end
            end
        end
            
        CompiledParticleSchnitzMap = [CompiledParticles{ChN}(:).Nucleus]; 
        if length(CompiledParticleSchnitzMap) > length(unique(CompiledParticleSchnitzMap))
            DoubledSchnitzes = [];
            for sc_idx = unique(CompiledParticleSchnitzMap)
                if length(find(CompiledParticleSchnitzMap == sc_idx)) > 1
                    DoubledSchnitzes = [DoubledSchnitzes, sc_idx];
                end
            end
            CPsToRemove = [];
            for sc_idx = DoubledSchnitzes
                CPlist = find(CompiledParticleSchnitzMap == sc_idx);
                ApprovedDoubledSchnitz = zeros(1,length(CPlist)); 
                NumFrames = zeros(1,length(CPlist)); 
                for j = 1:length(CPlist)
                    cp_idx = CPlist(j);
                    ApprovedDoubledSchnitz(j) = CompiledParticles{ChN}(cp_idx).Approved;
                    NumFrames(j) = length(CompiledParticles{ChN}(cp_idx).Frame);
                end
                if sum(ApprovedDoubledSchnitz > 0) == 1
                    for j = 1:length(CPlist)
                        if ApprovedDoubledSchnitz(j) < 1
                            CPsToRemove = [CPsToRemove, CPlist(j)];
                        end
                    end
                else
                    for j = 1:length(CPlist)
                        if NumFrames(j) < max(NumFrames)
                            CPsToRemove = [CPsToRemove, CPlist(j)];
                        end
                    end
                end
            end
            for i = 1:length(CPsToRemove);
                cp_idx = CPsToRemove(i);
                CompiledParticles{ChN}(cp_idx) = [];
                CPsToRemove = CPsToRemove-1;
            end
        end
    end
end


