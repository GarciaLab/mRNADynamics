function schnitzCellIndexes = schnitzCellsWithinFrames(Prefix, firstFrame,lastFrame)
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
% Author (contact): Gabriella Martini (martini@berkeley.edu)
% Created: 06/07/2020
% Last Updated: 06/07/2020

liveExperiment = LiveExperiment(Prefix);
[~,~,DropboxFolder,~,~]=...
    DetermineLocalFolders(Prefix);

PreProcPath = liveExperiment.userPreFolder;

DataFolder = [DropboxFolder, filesep, Prefix];
FilePrefix = [Prefix, '_'];
[FrameInfo, schnitzcells] =...
        loadCheckNuclearTrackingMats(DataFolder, PreProcPath, FilePrefix);
schnitzcellsSize = size(schnitzcells);

schnitzCellIndexes = [];
for i = 1:schnitzcellsSize(2)
    frameArraySize = size(schnitzcells(i).frames); % Might need to use Particles{CurrentChannel}(i).Frame
    inRange = 0;
    for k = 1:frameArraySize(2)
        currentFrame = schnitzcells(i).frames(k);
        if currentFrame >= firstFrame && currentFrame <= lastFrame
            inRange = 1;
        end
    end
    if inRange
        schnitzCellIndexes = [schnitzCellIndexes i];
    end
end

if isempty(schnitzCellIndexes)
    schnitzCellIndexes = -1; %This means no particles have been found in this range
end

end