function [Clusters] = sortSpotsByNuclei(prefix)

liveExperiment = LiveExperiment(prefix);

% I/O
dropboxFolder = liveExperiment.userResultsFolder;
OutputFolder = [DropboxFolder, filesep, Prefix];
experimentType = liveExperiment.experimentType;

inputChannels = liveExperiment.inputChannels;
spotChannels = liveExperiment.spotChannels;
clusterChannel = intersect(inputChannels, spotChannels);

if length(clusterChannel)~=1
    error('More than one input+spots (cluster) channel detected. This function is currently limited to processing one channel.')
end

% Iterate over all frames (b/c Spots structure is organized by frame)
for currFrame = 1:length(Spots{clusterChannel})
    
end