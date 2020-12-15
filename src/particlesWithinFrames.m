function particlesIndexes = particlesWithinFrames(Prefix, firstFrame,lastFrame, UseCompiledParticles, CurrentChannel)
% particlesWithinFrames(Prefix, Particles, firstFrame, lastFrame)
%
% DESCRIPTION
% Finds particles within the given frame range. For a particle to be in 
% range it needs to have at least one frame within the given range.
% 
% ARGUEMENTS
% Prefix: Prefix of the data set to analyze
% particles: An array of all particles in the movie.
% firstFrame: Start of the desired frame range  
% lastFrame: End of hte desired frame range
%
% OUTPUT
% This returns an array of the indexes of the particles found
% in the range or -1 if no particles were found.
%
% Author (contact): Emma Luu (emma_luu@berkeley.edu)
% Created: 06/05/2017
% Last Updated: 06/21/2017
%
% Documented by: Emma Luu (emma_luu@berkeley.edu)

if ~exist('UseCompiledParticles', 'var')
    UseCompiledParticles = false;
end

[~,~,DropboxFolder,~,~]=...
    DetermineLocalFolders(Prefix);
DataFolder=[DropboxFolder,filesep,Prefix];
if UseCompiledParticles
    load([DataFolder,filesep,'CompiledParticles.mat'], 'CompiledParticles');
    if isfield(CompiledParticles, 'CompiledParticles')
        CompiledParticles = CompiledParticles.CompiledParticles;
    end
    if iscell(CompiledParticles)
        Particles = CompiledParticles{CurrentChannel};
    else
        Particles = CompiledParticles;
    end
else
    load([DataFolder,filesep,'Particles.mat'])\
    particlesSize = size(Particles);
end


particlesIndexes = [];
for i = 1:particlesSize(2)
    frameArraySize = size(Particles(i).Frame); % Might need to use Particles{CurrentChannel}(i).Frame
    inRange = 0;
    for k = 1: frameArraySize(2)
        currentFrame = Particles(i).Frame(k);
        if currentFrame >= firstFrame && currentFrame <= lastFrame
            inRange = 1;
        end
    end
    if inRange
        particlesIndexes = [particlesIndexes i];
    end
end

if isempty(particlesIndexes)
    particlesIndexes = -1; %This means no particles have been found in this range
end

end