function Clusters = assignClusters2Nucleus(...
                            schnitzcells,Ellipses,Clusters,Spots,...
                            currFrame,pixelSizeXY,pixelSizeZ)
% DESCRIPTION
% Assign all spots to the closest nucleus. No restrictions on number of 
% spots that can be affiliated with a nucleus.
% This funciton does not track or connect the clusters.
% This is written assuming that you segemented the spots on nucleus-masked
% data, i.e. there is previous enforcement that the clusters must be inside
% a nucleus boundary. 
%
% ARGUMENTS
% schnitzcells:
% Ellipses:
% Spots:
% currFrame:
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
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 2022-05-09
% Last Updated: 2022-06-01

%% Get spot positions

% From here on, "new" means "current frame" for both the nuclei and
% particles.
[newSpotsX, newSpotsY, newSpotsZ, newSpotsXYZ] = ...
                            getSpotsXYZCentroids3DGauss(Spots(currFrame));

% if there are no spots in this frame, no need to continue
if isempty(newSpotsX)
    return
end
    
%% Get nucleus positions. 
% The order of both the nuclei and particle position is raw (coming from 
% Spots and Ellipses, respectively), not the one in the Particles and 
% Nuclei structures.

% Determine if we have z position info for nuclei
ellipsesColumnCount = size(Ellipses{currFrame},2);
if ellipsesColumnCount == 10
    % Using Ellipses 10th column as z coordinates for nuclei
    % This is only true if nuclei were tracked in 3D
    hasNuclearZCoordinates = true;
    ellipsesColumns = [1,2,10];
else
    % Typically, nuclei are tracked in 2D from projected His frames
    hasNuclearZCoordinates = false;
    ellipsesColumns = [1,2];
end

%Get the position of the nuclei in this frame
newEllipsesXYZ = Ellipses{currFrame}(:,ellipsesColumns);
newEllipsesX = newEllipsesXYZ(:,1);
newEllipsesY = newEllipsesXYZ(:,2);
if hasNuclearZCoordinates
    newEllipsesZ = newEllipsesXYZ(:,3);
end

%% Calculate distances between the spots and the ellipses in this frame
for currSpot=1:length(newSpotsX)

    if hasNuclearZCoordinates
        newSpotPos = [newSpotsX(currSpot)*pixelSizeXY, newSpotsY(currSpot)*pixelSizeXY, newSpotsZ(currSpot)*pixelSizeZ];
        newEllipsePos =[newEllipsesX*pixelSizeXY, newEllipsesY*pixelSizeXY, newEllipsesZ*pixelSizeZ];
    else
        newSpotPos = [newSpotsX(currSpot)*pixelSizeXY, newSpotsY(currSpot)*pixelSizeXY];
        newEllipsePos =[newEllipsesX*pixelSizeXY, newEllipsesY*pixelSizeXY];
    end

    distances(currSpot,:) = vecnorm(newSpotPos - newEllipsePos, 2, 2);

end

%% Find the closest ellipse to each spot

% minDistEllipseIndex is a row vector. The position in minDistIndex
% corresponds to each spot. The value at that position corresponds to
% the closest ellipse.
if size(distances,2)>1  % Check whether we had multiple ellipses
    [minDistValue, minDistEllipseIndex] = min(distances',[],1);
else  % If we only had one Ellipse, we need to be careful
    minDistValue = distances;
    minDistEllipseIndex = ones(size(distances));
end

%% Identify the nucleus lineage to which the closest ellipse belongs

% minDistSchnitz is now a row vector. The position in 
% minDistSchnitz corresponds to each spot. The value at that 
% position corresponds to the closest schnitzcell.
minDistSchnitz = zeros(size(minDistEllipseIndex));

for currMinIndex = 1:length(minDistEllipseIndex)
    currMinEllipse = minDistEllipseIndex(currMinIndex);
    % Otherwise, find the schnitzcell whose Ellipse in the 
    % currentFrame (encoded as "cellno" in the Ellipses struct) is 
    % closest to the Spot in minDistIndex we're currently checking.
    for currSchnitzID = 1:length(schnitzcells)
        schnitz2Check = schnitzcells(currSchnitzID);
        ellipse2Check = schnitz2Check.cellno(schnitz2Check.frames==currFrame);

        if ellipse2Check==currMinEllipse
            minDistSchnitz(currMinIndex) = currSchnitzID;
        end
    end
end

%% Sort spot data by nucleus and add to Cluster structure 

% Get the schnitzcells that already have Spots assigned to them in Clusters
assignedSchnitz = [Clusters.Nucleus];

% Find all the schnitzcells that were closest to spots in this frame
uniqueMinDistSchnitz = unique(minDistSchnitz);

% By nucleus, add all closest spots into the Cluster structure
for n = 1:numel(Clusters)%length(uniqueMinDistSchnitz)
    % Check if we have spots associated with this nucleus
    currSchnitzID = Clusters(n).Nucleus;
    schnitzMatchIndex = find(uniqueMinDistSchnitz==currSchnitzID);
    
    if isempty(schnitzMatchIndex)
        continue
    else
        % Set emptyFlag to false to indicate we have clusters at this frame
        Clusters(n).ClustersByFrame(currFrame).emptyFlag = 0;
        
        % Grab the spots we need to assign to this nucleus
        spotToClustersIndex = find(minDistSchnitz==currSchnitzID);
        %Get the spots' XYZ positions
        spotToClustersX = newSpotsX(spotToClustersIndex);
        spotToClustersY = newSpotsY(spotToClustersIndex);
        spotToClustersZ = newSpotsZ(spotToClustersIndex);

        % Now we can add all the data to the giant, ugly nested structure
        for s = 1:numel(spotToClustersIndex)
            Clusters(n).ClustersByFrame(currFrame).ClusterFits(s).SpotsIndex = spotToClustersIndex(s);
            Clusters(n).ClustersByFrame(currFrame).ClusterFits(s).xPos = spotToClustersX(s);
            Clusters(n).ClustersByFrame(currFrame).ClusterFits(s).yPos = spotToClustersY(s);
            Clusters(n).ClustersByFrame(currFrame).ClusterFits(s).zPos = spotToClustersZ(s);
        end 
    end
end