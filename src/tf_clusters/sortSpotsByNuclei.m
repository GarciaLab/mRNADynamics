function Clusters = sortSpotsByNuclei(prefix)
%
% DESCRIPTION
% This function assigns all spots detected in an input"channel (i.e. TF 
% clusters) to the closest nucleus lineage. It does NOT do any tracking of
% the clusters. It requires that nuclear segmentation and tracking be
% complete, but does not rely on any particle tracking code.
% It is currently limited to a single nuclear cycle (hard-coded to nc14).
%
% ARGUMENTS
% prefix: prefix of the data set to analyze
%
% OPTIONS
% NA
%
% OUTPUT
% Clusters: structure organized by Nucleus, listing all clusters
%           affiliated with each nucleus for each time frame. 
%           Data fields include:
%               - Nucleus: each nucleus lineage with clusters has its own
%                          element in the structure
%               - Frames: all frames containing clusters
%               - SpotsIndex: index of each cluster as found in the 
%                             Spots.mat data structure
%               - xPos, yPos, zPos: xyz position of each cluster in each
%                                   frame
% Clusters.mat: file storing the Clusters structure, saved to the
%               ResutlsFolder for the prefix
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 2022/05/09
% Last Updated:

%% I/O
liveExperiment = LiveExperiment(prefix);
preProcFolder = liveExperiment.preFolder;
dropboxFolder = liveExperiment.userResultsFolder;
outputFolder = [dropboxFolder, filesep, prefix];

pixelSize = liveExperiment.pixelSize_um;
pixelZSize = liveExperiment.zStep_um;

anaphaseFrames = liveExperiment.anaphaseFrames';
nc14 = anaphaseFrames(6);

inputChannels = liveExperiment.inputChannels;
spotChannels = liveExperiment.spotChannels;
nSpotCh = length(spotChannels);

clusterChannel = intersect(inputChannels, spotChannels);
nClusterCh = length(clusterChannel);
if nClusterCh~=1
    error('More than one input+spots (cluster) channel detected. This function is currently limited to processing one channel.')
end

%% Load nuclear segmentation/tracking, spots, and cluster data
[Ellipses, schnitzcells] = loadNucleiSegmentationAndLineages(outputFolder, prefix);
[Spots, ~] = loadSpotsAndCreateSpotFilter(dropboxFolder, prefix, nSpotCh);
Spots = Spots{1,clusterChannel};
% For now, just grab the one structure for the cluster channel to make
% indexing cleaner. 
% TODO: add a loop over all channels
Clusters = struct('Nucleus',{},'Frames',{},'SpotsIndex',{},...
                  'xPos',{},'yPos',{},'zPos',{});

%% Iterate over all frames (b/c Spots structure is organized by frame)
h = waitbar(0);
% Right now, we're limiting this to nc14 only
for currFrame = nc14:length(Spots)
    waitbar((currFrame-nc14)/(length(Spots)-nc14),h, ...
            sprintf('Sorting clusters from frame: %d of %d', currFrame, length(Spots)));
    Clusters = assignClusters2Nucleus(schnitzcells,Ellipses,Clusters,...
                                      Spots,currFrame,pixelSize,pixelZSize);
%     disp(['processing frame:', num2str(currFrame)]);
                                  
end
close(h)

%% Save output
save([outputFolder, filesep, 'Clusters.mat'], 'Clusters');
