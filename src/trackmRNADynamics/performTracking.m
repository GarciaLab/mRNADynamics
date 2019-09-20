function [Particles, SpotFilter] = performTracking(Particles, scurrentChannelitzcells, NCh, Spots, app, SpotFilter, PreProcPath, Prefix, UseHistone, ParticlesFig, SpotsChannel, NucleiFig, particlesAxes, nucAxes, Ellipses, PixelSize, SearchRadius, ExperimentType, FrameInfo, retrack, displayFigures)


NDigits = adjustIndexSizeAccordingToFrames(FrameInfo);


% Iterate over all channels

for Channel = 1:NCh
    
    % Iterate over all frames
    for CurrentFrame = 1:length(Spots{Channel})
        
        if isempty(app) && displayFigures
            figure(ParticlesFig)
            set(ParticlesFig, 'units', 'normalized', 'position', [0.01, .55, .33, .33]);
        end
        
        % Get the filter for this frame
        CurrentFrameFilter = logical(SpotFilter{Channel}(CurrentFrame, ~isnan(SpotFilter{Channel}(CurrentFrame, :))));
        
        xPos = displayParticlesFigure(app, particlesAxes, ParticlesFig, Spots, Channel, CurrentFrame, ...
            CurrentFrameFilter, PreProcPath, Prefix, SpotsChannel, FrameInfo, displayFigures);
        
        if UseHistone
            [Particles, SpotFilter] = trackParticlesBasedOnNuclei(PreProcPath, Prefix, CurrentFrame, NDigits, app, nucAxes, Ellipses, ...
                ExperimentType, Channel, scurrentChannelitzcells, Particles, Spots, SpotFilter, PixelSize, SearchRadius, retrack, displayFigures);
        else
            [Particles] = trackParticlesBasedOnProximity(Particles, Spots, xPos, SpotFilter, Channel, CurrentFrame, PixelSize, SearchRadius, retrack, displayFigures);
        end
        
    end
    
end

if isempty(app) && displayFigures
    close(ParticlesFig)
    if UseHistone
        close(NucleiFig)
    end
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
    
    %AR: add x, y and z positions to the Particles structure. This was
    %originally only added by addParticlePosition, but it it's more useful
    %early on in the pipeline.
    
    for particle=1:length(Particles{currentChannel})
        for frame=1:length(Particles{currentChannel}(particle).Frame)
            [x,y,z]=SpotsXYZ(Spots{currentChannel}(Particles{currentChannel}(particle).Frame(frame)));
            if ~isempty(x)
                Particles{currentChannel}(particle).xPos(frame)=x(Particles{currentChannel}(particle).Index(frame));
                Particles{currentChannel}(particle).yPos(frame)=y(Particles{currentChannel}(particle).Index(frame));
                Particles{currentChannel}(particle).zPos(frame)=z(Particles{currentChannel}(particle).Index(frame));
            end
        end
    end
    
end

% If we only have one channel, then convert SpotFilter and Particles to a standard structure.
if NCh == 1
    SpotFilter = SpotFilter{1};
    Particles = Particles{1};
end

end


% See how  many frames we have and adjust the index size of the files to
% load accordingly
function NDigits = adjustIndexSizeAccordingToFrames(FrameInfo)

if length(FrameInfo) < 1E3
    NDigits = 3;
elseif length(FrameInfo) < 1E4
    NDigits = 4;
else
    error('No more than 10,000 frames supported. Change this in the code')
end

end
