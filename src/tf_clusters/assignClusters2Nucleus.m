function [Clusters] = assignClusters2Nucleus(...
                            schnitzcells,Ellipses,Particles,Spots,...
                            currFrame,pixelSize, pixelZSize)
% DESCRIPTION
% Assign all spots to the closest nucleus. No restrictions on number of 
% spots that can be affiliated with a nucleus.
% This funciton does not track or connect the clusters.
% This is written assuming that you segemented the spots on nucleus-masked
% data, i.e. there is previous enforcement that the clusters must be inside
% a nucleus boundary. 
%
% ARGUMENTS
% schnitzcells
% Ellipses
% Particles
% Spots
% currFrame
%
% OPTIONS
% NA
%
% OUTPUT
% Clusters.mat: structure organized by Nucleus, listing all spots
%               affiliated with each nucleus for each time frame
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 05/09/2022
% Last Updated:


%% Get the particles and nuclei positions. 
% The order of both the nuclei and particle position is raw (coming from 
% Spots and Ellipses, respectively), not the one in the Particles and 
% Nuclei structures.

% From here on, "curr" means "current frame" for both the nuclei and
% particles.
[currSpotsX,currSpotsY,currSpotsZ] = getSpotsXYZ(Spots(currFrame));

if ~isempty(currSpotsX)

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
    currEllipsesXYZ = Ellipses{currFrame}(:,ellipsesColumns);
    currEllipsesX = currEllipsesXYZ(:,1);
    currEllipsesY = currEllipsesXYZ(:,2);
    if hasNuclearZCoordinates
        currEllipsesZ = currEllipsesXYZ(:,3);
    end

    %Calculate the distances between the spots in this frame and the
    %ellipses
    for jSchnitz=1:length(currSpotsX)

        if hasNuclearZCoordinates
            currSpotPos = [currSpotsX(jSchnitz)*pixelSize, currSpotsY(jSchnitz)*pixelSize, currSpotsZ(jSchnitz)*pixelZSize];
            currEllipsePos =[currEllipsesX*pixelSize, currEllipsesY*pixelSize, currEllipsesZ*pixelZSize];
        else
            currSpotPos = [currSpotsX(jSchnitz)*pixelSize, currSpotsY(jSchnitz)*pixelSize];
            currEllipsePos =[currEllipsesX*pixelSize, currEllipsesY*pixelSize];
        end
        
        distances(jSchnitz,:) = vecnorm(currSpotPos - currEllipsePos, 2, 2);

    end
    
    % minDistIndex is a row vector. The position in minDistIndex
    % corresponds to each spot. The value at that position corresponds to
    % the closest ellipse.
    if size(distances,2)>1  % Check whether we had multiple ellipses
        [minDistValue,minDistIndex] = min(distances,[],1);
    else  % If we only had one Ellipse, we need to be careful
        minDistValue = distances;
        minDistIndex = ones(size(distances));
    end
    
    
    % Find the schnitzcell corresponding to minDistIndices in the 
    % schnitzcell structure.
    %
    % minDistSchnitz is now a row vector. The position in 
    % minDistSchnitz corresponds to each spot. The value at that 
    % position corresponds to the closest schnitzcell.
    minDistSchnitz = zeros(size(minDistIndex));
    
    for iDistIndex = 1:length(minDistIndex)
        distIndex2Check = minDistIndex(iDistIndex);
        
        % If retrack, set the minDistSchnitz of the already assigned spots
        % to infinity
        if minDistIndex(iDistIndex)==inf
            minDistSchnitz(iDistIndex) = inf;
        else
            % Otherwise, find the schnitzcell whose Ellipse in the 
            % currentFrame (encoded as "cellno" in the Ellipses struct) is 
            % closest to the Spot in minDistIndex we're currently checking.
            for jSchnitz = 1:length(schnitzcells)
                schnitz2Check = schnitzcells(jSchnitz);
                ellipse2Check = schnitz2Check.cellno(schnitz2Check.frames==currFrame);
                
                if ellipse2Check==distIndex2Check
                    minDistSchnitz(iDistIndex) = jSchnitz;
                end
            end
        end
    end

   
    %Find the schnitz assigned to each already existing particle
    AssignedSchnitz=[];
    for i=1:length(Particles)
        if ~isempty(Particles(i).Nucleus)
            AssignedSchnitz(i)=Particles(i).Nucleus;
        else
            AssignedSchnitz(i)=0;
        end
    end

    
    %Find the schnitz that were closest to the spots found in this frame.
    UniqueMinIndexSchnitz=unique(minDistSchnitz);
    for iDistIndex=1:length(UniqueMinIndexSchnitz)
      
        
       %Find the particles that are assigned to this schnitz
        ParticleToAssign=find(AssignedSchnitz==UniqueMinIndexSchnitz(iDistIndex));

        %Are there any particles in previous frames that are assigned to
        %this schnitz? If not, we move on and define the current spots as
        %new particles.
        if ~isempty(ParticleToAssign)
                        
            %Find the spots I need to locate to this schnitz
            SpotToParticleIndices=find(minDistSchnitz==UniqueMinIndexSchnitz(iDistIndex));
            %Get the spots' XY positions
            SpotToParticleX=currSpotsX(SpotToParticleIndices);
            SpotToParticleY=currSpotsY(SpotToParticleIndices);

            %Positions of the particles. Note that I'm not enforcing
            %that the Particle have a spot in the previous frame. I'll
            %have to see whether this leads ot any issues.
            PreviousParticlesX=[];
            PreviousParticlesY=[];
            for jSchnitz=1:length(ParticleToAssign)
                [PreviousSpotsX,PreviousSpotsY]=...
                    getSpotsXYZ(Spots(Particles(ParticleToAssign(jSchnitz)).Frame(end)));
                PreviousParticlesX(jSchnitz)=...
                    PreviousSpotsX(Particles(ParticleToAssign(jSchnitz)).Index(end));
                PreviousParticlesY(jSchnitz)=...
                    PreviousSpotsY(Particles(ParticleToAssign(jSchnitz)).Index(end));
            end

            %Calculate the distance of the Spots to the Particles
            %within this schnitz.

            

            
            %Get the spots that are the nearest neighbors to
            %each the previous particles
            [NewSpotToAssign,DistancesToNewSpots] =...
                knnsearch([SpotToParticleX;SpotToParticleY]'*pixelSize,...
                [PreviousParticlesX;PreviousParticlesY]'*pixelSize,...
                'K',SpotsPerNucleus);
            
            %Each row in NewSpotToAssign corresponds to the previous
            %particles. The index within that row tells us which new spot
            %it is closest to. Subsequent columns go beyond the nearest
            %neighbor, and tell us about higher-order neighbors.
            
%             %If the number of spots less than SpotsPerNucleus, then the
%             %code doesn't generate the 2nd and higher order nearest
%             %neighbors. In that case, pad the results from knnsearch.
%             if size(DistancesToNewSpots,2)<SpotsPerNucleus
%                 DistancesToNewSpots=padarray(DistancesToNewSpots,[0,SpotsPerNucleus-size(DistancesToNewSpots,2)],'post');
%                 NewSpotToAssign=padarray(NewSpotToAssign,[0,SpotsPerNucleus-size(NewSpotToAssign,2)],'post');
%                 DistancesToNewSpots(DistancesToNewSpots==0)=nan;
%                 NewSpotToAssign(NewSpotToAssign==0)=nan;
%             end
            
            
            %Sometimes, two or more new spots are nearest to the same particle. We want to make
            %sure that the farther one gets assigned to the next most appropriate
            %particle.
            %I'll go though each unique set of indices in each column of
            %NewSpotToAssign (which correspond to the new spots).
            for jSchnitz=1:size(NewSpotToAssign,2)
                %Find the indices that are repeated in this column. This means that we
                %have multiple new spots assigned to the previous particle.
                UniqueNewSpotIndices=unique(NewSpotToAssign(:,jSchnitz));
                for k=1:length(UniqueNewSpotIndices)
                    %Find the positions in this column that correspond to the multiply
                    %assigned new spots
                    Positions=find(NewSpotToAssign(:,jSchnitz)==UniqueNewSpotIndices(k));
                    %Find the minimum distance, and set the other ones to infinity
                    [~,MinIndexPositions]=min(DistancesToNewSpots(Positions,jSchnitz));
                    DistancesToNewSpots(Positions(~ismember(Positions,Positions(MinIndexPositions))),jSchnitz)=inf;
                    %Also, set the distance to this spot with respect to
                    %the other particles to infinity.
                    DistancesToNewSpots(Positions(ismember(Positions,Positions(MinIndexPositions))),jSchnitz+1:end)=inf;
                end
            end
            %Now, we can remake the assignment list and the distances
            [DistancesToNewSpots,NewSpotToAssignIndices]=min(DistancesToNewSpots,[],2);
            for jSchnitz=1:length(NewSpotToAssignIndices)
                NewSpotToAssignTemp(jSchnitz)=NewSpotToAssign(jSchnitz,NewSpotToAssignIndices(jSchnitz));
            end
            NewSpotToAssign=NewSpotToAssignTemp';
            

            %Assign the new spots to their corresponding particles
            for jSchnitz=1:length(ParticleToAssign)
                if ~isinf(DistancesToNewSpots(jSchnitz))
                    SpotIndexToCopy=SpotToParticleIndices(NewSpotToAssign(jSchnitz));
                    %Finally, copy the information onto this particle.
                    if retrack == 0 || Particles(ParticleToAssign(jSchnitz)).Approved < 1
                        Particles(ParticleToAssign(jSchnitz)).Index(end+1)=SpotIndexToCopy;
                        Particles(ParticleToAssign(jSchnitz)).Frame(end+1)=currFrame;
                    end
                    %Remove this spot from the pool so that it doesn't get
                    %assigned to a new particle at the end of the code
                    newParticlesFlag(SpotIndexToCopy)=0;
                end
            end
            
        end
    end


    %Assign the particles that correspond to a new nucleus
    NewParticlesIndices=find(newParticlesFlag==1);

    %Indices of the particles that won't be assigned to a nucleus and
    %therefor need to be disapproved by flagging them in SpotFilter
    IndexToMove=[];

   
    for iDistIndex=1:length(NewParticlesIndices)

        %Recalculate the assigned nuclei in this frame. Note that this is different
        %from AssignedSchnitz in the sense that it asks how many Particles exist in
        %this frame associated with a given Schnitz. 
        %This is useful in case two new particles are closed to a given new nucleus. Right now it
        %will assign the first particle found. I might have to change
        %this if it becomes too annoying.
        AssignedSchnitzCurrentFrame=[];         %These keeps track of which nuclei have
                                    %already been assigned to particles
        for jSchnitz=1:length(Particles)
            if sum(Particles(jSchnitz).Frame==currFrame)
                AssignedSchnitzCurrentFrame=[AssignedSchnitzCurrentFrame,Particles(jSchnitz).Nucleus];
            end
        end


        %Make sure the total number of particles assigned to this schnitz
        %is not higher than SpotsPerNucleus
        if sum(AssignedSchnitzCurrentFrame==minDistSchnitz(NewParticlesIndices(iDistIndex)))<SpotsPerNucleus
            Particles(end+1).Frame=currFrame;
            Particles(end).Index=NewParticlesIndices(iDistIndex);
            Particles(end).Nucleus=minDistSchnitz(NewParticlesIndices(iDistIndex));
            
            assert(Particles(end).Nucleus <= length(schnitzcells));
            
            if retrack
                Particles(end).Approved=0;
            end
        else
            %If the spot cannot be assigned to a particle and schnitz, then we'll
            %take the spot out using SpotFilter.
            IndexToMove=[IndexToMove,iDistIndex];
        end
    end

    if ~isempty(IndexToMove)     
        SpotFilter(currFrame,IndexToMove)=0;
    end


end