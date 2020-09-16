function [Particles, SpotFilter] =...
    performTracking(Particles, schnitzcells, NCh, Spots, app,...
    SpotFilter, PreProcPath, Prefix, UseHistone, ParticlesFig,...
    SpotsChannel, NucleiFig, particlesAxes, nucAxes,...
    Ellipses, PixelSize, SearchRadius, ExperimentType,...
    FrameInfo, retrack, displayFigures, thisExperiment)

NDigits = thisExperiment.nDigits;


% Iterate over all channels
h = waitbar(0, 'Performing tracking ...');
for Channel = 1:NCh
    
    % Iterate over all frames

    for CurrentFrame = 1:length(Spots{Channel})
        waitbar(CurrentFrame/length(Spots{Channel}),h,'Performing tracking ...');
        if isempty(app) && displayFigures
            figure(ParticlesFig)
            set(ParticlesFig, 'units', 'normalized', 'position', [0.01, .55, .33, .33]);
        end
        
        % Get the filter for this frame
        CurrentFrameFilter = logical(SpotFilter{Channel}(CurrentFrame,...
            ~isnan(SpotFilter{Channel}(CurrentFrame, :))));
        
        xPos = displayParticlesFigure(app, particlesAxes,...
            ParticlesFig, Spots, Channel, CurrentFrame, ...
            CurrentFrameFilter, PreProcPath, Prefix, SpotsChannel, FrameInfo, displayFigures);
        
        if UseHistone
            [Particles, SpotFilter] = trackParticlesBasedOnNuclei(PreProcPath, Prefix,...
                CurrentFrame, NDigits, app, nucAxes, Ellipses, ...
                ExperimentType, Channel, schnitzcells, Particles, Spots,...
                SpotFilter, PixelSize, SearchRadius, retrack, displayFigures);
        else
            [Particles] = trackParticlesBasedOnProximity(Particles, Spots,...
                xPos, SpotFilter, Channel, CurrentFrame, PixelSize, SearchRadius,...
                retrack, displayFigures);
        end
        
    end
    
end

if isempty(app) && displayFigures
    close(ParticlesFig)
    if UseHistone, close(NucleiFig); end
end


for currentChannel = 1:NCh
    
    if ~isfield(Particles{currentChannel}, 'FrameApproved')
        
        for particle = 1:length(Particles{currentChannel})
            Particles{currentChannel}(particle).FrameApproved = true(size(Particles{currentChannel}(particle).Frame));
        end
        
    else
        
        for particle = 1:length(Particles{currentChannel})
            
            if isempty(Particles{currentChannel}(particle).FrameApproved)
                Particles{currentChannel}(particle).FrameApproved = true(size(Particles{currentChannel}(particle).Frame));
            end
            
        end
        
    end
    
    
    Particles = addPositionsToParticles(Particles, Spots, currentChannel);
    
end

% If we only have one channel, then convert SpotFilter and Particles to a standard structure.
if NCh == 1
    SpotFilter = SpotFilter{1};
    Particles = Particles{1};
end


close(h)
end