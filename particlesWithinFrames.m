function particlesIndexes = particlesWithinFrames(Prefix, Particles, firstFrame,lastFrame)
% Finds particles within a given frame range and returns an array of their
% indexes
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);
DataFolder=[DropboxFolder,filesep,Prefix];
load([DataFolder,filesep,'Particles.mat']);
particlesSize = size(Particles);

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
end