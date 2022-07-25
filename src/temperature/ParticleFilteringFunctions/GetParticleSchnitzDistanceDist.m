function ParticleSchnitzDistances = GetParticleSchnitzDistanceDist(schnitzcells, CompiledParticles, xDim, yDim)
%% Calculates all particle schnitz distances for schnitz/particle pairs that are well within the embryo boundaries
if iscell(CompiledParticles)
    CompiledParticles = CompiledParticles{1};
end
ParticleSchnitzDistances = [];
for p = 1:length(CompiledParticles)
    CurrentParticle = CompiledParticles(p);
    CurrentSchnitz = schnitzcells(CurrentParticle.schnitz);
    if (CurrentSchnitz.Approved == 1) & (CurrentParticle.Approved == 1)
        schnitzframes = CurrentSchnitz.frames.';
        particleframes = CurrentParticle.Frame;
        for sc_index = 1:length(schnitzframes)
            if ((CurrentSchnitz.FrameApproved(sc_index)==1) && ...
                    (CurrentSchnitz.cenx(sc_index) > xDim/8) && (CurrentSchnitz.cenx(sc_index) < (xDim-xDim/8)) && ...
                    (CurrentSchnitz.ceny(sc_index) > yDim/8) && (CurrentSchnitz.ceny(sc_index) < (yDim-yDim/8)))
                if ~isempty(find( particleframes == schnitzframes(sc_index)))
                    p_index = find( particleframes == schnitzframes(sc_index));
                    if CurrentParticle.FrameApproved(p_index) == 1
                        ParticleSchnitzDistances(end+1) = sqrt((CurrentSchnitz.cenx(sc_index) -CurrentParticle.xPos(p_index))^2 + ...
                            (CurrentSchnitz.ceny(sc_index) -CurrentParticle.yPos(p_index))^2);
                    end
                end
                
            end
        end
    end
end