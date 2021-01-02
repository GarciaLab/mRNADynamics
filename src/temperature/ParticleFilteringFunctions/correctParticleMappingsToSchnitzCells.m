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
    try
        CompiledParticles = getCompiledParticles(liveExperiment);
    end
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
            end
        end
    end
end


