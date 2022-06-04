function ms2ClusterDistances = pairwiseDistMS2Clusters(prefix, varargin)
%
% DESCRIPTION
% This function takes segmented MS2 spots and TF clusters (using the
% Spots.m structure) and calculates the pairwise distances between the 
% MS2 spot in a nucleus and all clusters detected in the same nucleus.
%
% Note: This function currently only works for 1spot data
% Note: This function defaults to only processing nc14, unless otherwise
%       specified by the input options
%
%
% INPUT ARGUMENTS
% Prefix: prefix for this experiment
% 
% OPTIONS
% 'nc', #: 
%
% OUTPUT
% ms2ClusterDistances: A structure containing the pairwise distances  
%                      between an MS2 spot and all TF clusters detected,
%                      for each nucleus and time frame.
%
% ms2ClusterDistances.mat: MATLAB file containing the ms2ClusterDistances
%                          structure, saved to a subdirectory, entitled
%                          'clusterAnalysis', within your Results folder
%                          (aka, your 'DropboxFolder')
%
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 2020/08/12
% Last Updated: 2022/05/11
%

%% Parse input arguments
liveExperiment = LiveExperiment(prefix);

ncToAnalyze = 14;
for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'nc')
        ncToAnalyze = varargin{i+1};
    end
end

% Set nc start time
anaphaseFrames = liveExperiment.anaphaseFrames';
nc13 = anaphaseFrames(5);
nc14 = anaphaseFrames(6);

if ncToAnalyze==13
    ncStart = nc13;
else
    ncStart = nc14;
end
    

%% Get info from LiveExperiment
nFrames = liveExperiment.nFrames;
% pixel sizes in um for accurate distance calculations
pixelSizeXY = liveExperiment.pixelSize_um;
pixelSizeZ = liveExperiment.zStep_um;

% channel info
inputChannels = liveExperiment.inputChannels;
spotChannels = liveExperiment.spotChannels;
nSpotCh = numel(spotChannels);

clusterChannel = intersect(inputChannels, spotChannels);
nClusterCh = length(clusterChannel);

ms2Channel = setxor(inputChannels, spotChannels);
nMs2Ch = length(ms2Channel);
% enforce 1 cluster + 1 MS2 channel
if nClusterCh~=1
    error('More than one input+spots (cluster) channel detected. This function is currently limited to processing one cluster input channel.')
elseif nMs2Ch~=1
    error('More than one MS2 spot channel detected. This function is currently limited to processing one transcriptional output channel.')
end

% I/O
dropboxFolder = liveExperiment.userResultsFolder;
resultsFolder = [dropboxFolder, filesep, prefix];
writeFolder = [resultsFolder filesep 'cluster_analysis'];
if ~exist(writeFolder,'dir')
    mkdir(writeFolder)
end

%% Load cluster & particle data
[Spots, ~] = loadSpotsAndCreateSpotFilter(dropboxFolder, prefix, nSpotCh);
Spots = Spots{1,ms2Channel}; %we only need the MS2 channel here

Clusters = loadClusters(resultsFolder);

[Particles, ~] = getParticles(liveExperiment);
Particles = Particles{1,ms2Channel};


%% Quality control
% QC for MS2 traces
% minTraceLen = ceil((nFrames - ncStart)/2);
minTraceLen = 15; % 3 min
Particles = filterParticlesForClusterAnalysis(Particles, ncStart, minTraceLen);

% QC for nuclear traces
[~, schnitzcells] = loadNucleiSegmentationAndLineages(resultsFolder, prefix);
Clusters =  filterClustersForClusterAnalysis(Clusters, schnitzcells);

%% Get only nuclei that have both an MS2 trace & clusters detected
nucleiWithMS2 = [Particles.Nucleus];
if length(nucleiWithMS2)~= length(unique(nucleiWithMS2))
    error('Non-unique nuclei in the MS2 Particles structure - this function is not designed to handle this scenario.')
end
nucleiWithClusters = [Clusters.Nucleus];
nucleiWithBoth = intersect(nucleiWithMS2, nucleiWithClusters);
numNucleiWithBoth = numel(nucleiWithBoth);

%% Loop over all nuclei with both MS2 spots & TF clusters
ms2ClusterDistances(numNucleiWithBoth).nucleusID = [];  %pre-allocated for speed
for n = 1:numNucleiWithBoth
    currNucID = nucleiWithBoth(n);
    ms2ClusterDistances(n).nucleusID = currNucID;
    currClusters = Clusters(find(nucleiWithClusters == currNucID)).ClustersByFrame; %#ok<FNDSB>
    currParticle = Particles(find(nucleiWithMS2 == currNucID)); %#ok<FNDSB>
    
    % Go through each frame
    for f = 1:nFrames
        ms2ClusterDistances(n).frames(f).frameID = f;
        % initializing all fields in case there are no clusters/MS2 spots
        % in this frame
        ms2ClusterDistances(n).frames(f).xyzCoordClusters = [];
        ms2ClusterDistances(n).frames(f).xyzCoordMS2 = [];
        ms2ClusterDistances(n).frames(f).ms2ClusterDist = [];
        xyzPosMS2 = [];
        xyzPosClusters = [];
        
        %% Figure out whether this nucleus has both an MS2 spot and 
        % cluster(s) in this frame (less likely to have an MS2 spot)
        frameHasClusters = ~currClusters(f).emptyFlag; % flip so true if clusters are in frame
        ms2FrameIndex = find(currParticle.Frame == f);
        
        if ~isempty(ms2FrameIndex) && frameHasClusters
            %% Get all the xyz positions for the clusters and MS2 spot
            % putting the xyz coords in a structure that will be easy 
            % to pass directly to pdist2 later on in this script
            
            % The xyz positions stored in Particles are from 2D fits, not
            % the 3D Gaussian fits. Have to go back to Spots to get those
            spotsIndex = currParticle.Index(ms2FrameIndex);
            ms2Pos3D = Spots(f).Fits(spotsIndex).GaussPos3D;
            
            % Positional data is in pixels, need to convert to nm to 
            % correctly calculate distances.
            xyzPosMS2 = [ms2Pos3D(1)*pixelSizeXY,...
                         ms2Pos3D(2)*pixelSizeXY,...
                         ms2Pos3D(3)*pixelSizeZ];
            xyzPosClusters(:,1) = [currClusters(f).ClusterFits.xPos]*pixelSizeXY;
            xyzPosClusters(:,2) = [currClusters(f).ClusterFits.yPos]*pixelSizeXY;
            xyzPosClusters(:,3) = [currClusters(f).ClusterFits.zPos]*pixelSizeZ;
            
%             % This is the wrong way to do it, b/c assumes isometric xyz
%             % pixel size, when actually z is much larger than xy
%             xyzPosMS2 = [ms2Pos3D(1),...
%                          ms2Pos3D(2),...
%                          ms2Pos3D(3)];
%             xyzPosClusters(:,1) = [currClusters(f).ClusterFits.xPos];
%             xyzPosClusters(:,2) = [currClusters(f).ClusterFits.yPos];
%             xyzPosClusters(:,3) = [currClusters(f).ClusterFits.zPos];
            
            % Remove duplicate spot detections to minimize redundancy
            % (not sure why these aren't getting filtered out in
            % segmentSpots ... MT: TODO)
            uniqueXYZPosClusters = unique(xyzPosClusters,'stable','rows');
            
            % Store the coordinates by nucleus, frame
            ms2ClusterDistances(n).frames(f).xyzCoordClusters = uniqueXYZPosClusters;
            ms2ClusterDistances(n).frames(f).xyzCoordMS2 = xyzPosMS2;

            ms2Coord = xyzPosMS2;
            clusterCoords = uniqueXYZPosClusters;

            %% Calculate Pairwise distances
            if ~isempty(ms2Coord) && ~isempty(clusterCoords)
                ms2ClusterDist = pdist2(ms2Coord, clusterCoords);
                ms2ClusterDistances(n).frames(f).ms2ClusterDist = ms2ClusterDist;
            else
                ms2ClusterDistances(n).frames(f).ms2ClusterDist = [];   % make sure it exists even if there's nothing there
            end
        end
    end
end

%% Save
save([writeFolder filesep 'ms2ClusterDistances_px.mat'],'ms2ClusterDistances');