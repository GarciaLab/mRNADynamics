function [Particles,SpotFilter]=AssignParticle2Nucleus(...
    schnitzcells,Ellipses,Particles,Spots,SpotFilter,...
    CurrentFrame,PixelSize,SpotsPerNucleus, retrack)

%Find the N closest particles to each nucleus. N is determined by
%ExperimentType (1spot, 2spots). If N>1, then the code also uses the spot
%position for tracking.


%Get the particles and nuclei positions. The order of both the nuclei and
%particle position is raw (coming from Spots and Ellipses, respectively),
%not the one in the Particles and Nuclei structures.

%From here on, "New" means "CurrentFrame" for both the nuclei and
%particles.
[NewSpotsX,NewSpotsY]=SpotsXYZ(Spots(CurrentFrame));

if ~isempty(NewSpotsX)

    %Get the position of the nuclei in this frame
    NewEllipsesXY=Ellipses{CurrentFrame}(:,[1,2]);
    NewEllipsesX=NewEllipsesXY(:,1);
    NewEllipsesY=NewEllipsesXY(:,2);

    %This is to keep track of already-assigned particles
    NewParticlesFlag=ones(size(NewSpotsX));


    %Calculate the distances between the spots on this frame and the
    %ellipses

    for j=1:length(NewSpotsX)

        newSpotPos = [NewSpotsX(j), NewSpotsY(j)]*PixelSize;
        newEllipsePos =[NewEllipsesX, NewEllipsesY]*PixelSize;
        
        Distance(j,:) = vecnorm(newSpotPos - newEllipsePos, 2, 2);

    end
    
        
    %If retrack make the distance for the particles that have been
    %approved infinite.
    if retrack
        for i=1:length(Particles)
            %Find which approved particles are in this frame
            if sum(Particles(i).Frame==CurrentFrame)&(Particles(i).Approved~=0)
                %Make the distance infinite
                Distance(Particles(i).Index(find(Particles(i).Frame==CurrentFrame)),:)=inf;
                NewParticlesFlag(Particles(i).Index(find(Particles(i).Frame==CurrentFrame)))=-1;
            end
        end
    end
    
    %MinIndex is a row vector. The position in MinIndex
    %corresponds to each spot. The value at that
    %position corresponds to the closest ellipse.
    if size(Distance,2)>1  %Check whether we had multiple ellipses
        [MinValues,MinIndex]=min(Distance');
    else %If we only had one Ellipse, we need to be careful
        MinValues=Distance;
        MinIndex=ones(size(Distance));
    end
               
    %If retrack, set the distance of the already assigned spots to
    %infinity
    if retrack
        MinIndex(MinValues==inf)=inf;
    end
    
    %Find the schnitz corresponding to MinIndex in the schnitzcell structure.
    %MinIndexNuclei is now a row vector. The position in MinIndexNuclei
    %corresponds to each spot. The value at that position corresponds
    %to the closest schnitz.
    MinIndexSchnitz=zeros(size(MinIndex));
    for i=1:length(MinIndex)
        for j=1:length(schnitzcells)
            try

                sj = schnitzcells(j);

                if sj.cellno(sj.frames==CurrentFrame) == MinIndex(i)
                    MinIndexSchnitz(i)=j;
                end
                
            catch
                 %AR 3/31/2019- sometimes this errors and I couldn't discover
                %why. 
            end
        end
        if MinIndex(i)==inf
            MinIndexSchnitz(i)=inf;
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
    UniqueMinIndexSchnitz=unique(MinIndexSchnitz);
    for i=1:length(UniqueMinIndexSchnitz)
      
        
       %Find the particles that are assigned to this schnitz
        ParticleToAssign=find(AssignedSchnitz==UniqueMinIndexSchnitz(i));

        %Are there any particles in previous frames that are assigned to
        %this schnitz? If not, we move on and define the current spots as
        %new particles.
        if ~isempty(ParticleToAssign)
                        
            %Find the spots I need to locate to this schnitz
            SpotToParticleIndices=find(MinIndexSchnitz==UniqueMinIndexSchnitz(i));
            %Get the spots' XY positions
            SpotToParticleX=NewSpotsX(SpotToParticleIndices);
            SpotToParticleY=NewSpotsY(SpotToParticleIndices);

            %Positions of the particles. Note that I'm not enforcing
            %that the Particle have a spot in the previous frame. I'll
            %have to see whether this leads ot any issues.
            PreviousParticlesX=[];
            PreviousParticlesY=[];
            for j=1:length(ParticleToAssign)
                [PreviousSpotsX,PreviousSpotsY]=...
                    SpotsXYZ(Spots(Particles(ParticleToAssign(j)).Frame(end)));
                PreviousParticlesX(j)=...
                    PreviousSpotsX(Particles(ParticleToAssign(j)).Index(end));
                PreviousParticlesY(j)=...
                    PreviousSpotsY(Particles(ParticleToAssign(j)).Index(end));
            end

            %Calculate the distance of the Spots to the Particles
            %within this schnitz.

            

            
            %Get the spots that are the nearest neighbors to
            %each the previous particles
            [NewSpotToAssign,DistancesToNewSpots] =...
                knnsearch([SpotToParticleX;SpotToParticleY]'*PixelSize,...
                [PreviousParticlesX;PreviousParticlesY]'*PixelSize,...
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
            for j=1:size(NewSpotToAssign,2)
                %Find the indices that are repeated in this column. This means that we
                %have multiple new spots assigned to the previous particle.
                UniqueNewSpotIndices=unique(NewSpotToAssign(:,j));
                for k=1:length(UniqueNewSpotIndices)
                    %Find the positions in this column that correspond to the multiply
                    %assigned new spots
                    Positions=find(NewSpotToAssign(:,j)==UniqueNewSpotIndices(k));
                    %Find the minimum distance, and set the other ones to infinity
                    [~,MinIndexPositions]=min(DistancesToNewSpots(Positions,j));
                    DistancesToNewSpots(Positions(~ismember(Positions,Positions(MinIndexPositions))),j)=inf;
                    %Also, set the distance to this spot with respect to
                    %the other particles to infinity.
                    DistancesToNewSpots(Positions(ismember(Positions,Positions(MinIndexPositions))),j+1:end)=inf;
                end
            end
            %Now, we can remake the assignment list and the distances
            [DistancesToNewSpots,NewSpotToAssignIndices]=min(DistancesToNewSpots,[],2);
            for j=1:length(NewSpotToAssignIndices)
                NewSpotToAssignTemp(j)=NewSpotToAssign(j,NewSpotToAssignIndices(j));
            end
            NewSpotToAssign=NewSpotToAssignTemp';
            
%             %Old Version:
%             %Assign the new spots to the previous particles
%             for j=1:length(ParticleToAssign)
%                 if ~isinf(DistancesToNewSpots(j))
%                     SpotIndexToCopy=SpotToParticleIndices(NewSpotToAssign(j));
%                     %Finally, copy the information onto this particle.
%                     Particles(ParticleToAssign(j)).Index(end+1)=SpotIndexToCopy;
%                     Particles(ParticleToAssign(j)).Frame(end+1)=CurrentFrame;
%                     %Remove this spot from the pool so that it doesn't get
%                     %assigned to a new particle at the end of the code
%                     NewParticlesFlag(SpotIndexToCopy)=0;
%                 end
%             end
            

            %Assign the new spots to their corresponding particles
            for j=1:length(ParticleToAssign)
                if ~isinf(DistancesToNewSpots(j))
                    SpotIndexToCopy=SpotToParticleIndices(NewSpotToAssign(j));
                    %Finally, copy the information onto this particle.
                    if retrack == 0 || Particles(ParticleToAssign(j)).Approved < 1
                        Particles(ParticleToAssign(j)).Index(end+1)=SpotIndexToCopy;
                        Particles(ParticleToAssign(j)).Frame(end+1)=CurrentFrame;
                    end
                    %Remove this spot from the pool so that it doesn't get
                    %assigned to a new particle at the end of the code
                    NewParticlesFlag(SpotIndexToCopy)=0;
                end
            end
            
            

%             [MinValuesParticles,MinIndexParticles]=SortDistances([SpotToParticleX;SpotToParticleY]',...
%                ,PixelSize)

%             
%             %Calculate the distances between the spots on this frame and the
%             %ellipses
%             clear DistanceParticles
% 
%             for j=1:length(SpotToParticleX)
%                 DistanceParticles(j,:)=sqrt((SpotToParticleX(j)*PixelSize-...
%                     PreviousParticlesX*PixelSize).^2+...
%                     (SpotToParticleY(j)*PixelSize-...
%                     PreviousParticlesY*PixelSize).^2);
%             end
%             %MinIndexParticles is a row vector. The position in MinIndex
%             %corresponds to each spot. The value at that
%             %position corresponds to the closest particle. Note that this
%             %is not true if we have only one particle
%             [MinValuesParticles,MinIndexParticles]=min(DistanceParticles');
%             
%             
%             %Sometimes the two or more new Spots are closer to the same
%             %previous Particle. In that case, only keep the closest
%             %Spot-Particle pair and set the distance of the other spots to
%             %the same particle to infinity.
%             UniqueMinIndexParticles=unique(MinIndexParticles);
%             
%                         
%             %Assign the spots to their corresponding particles. I need to
%             %consider the case where there was only one previous particle.
%             if length(ParticleToAssign)>1
%                 %If we have more than one previous particle, I can make use
%                 %of the matrix Distance, of MinValuesParticles, and
%                 %MinIndexParticles.
%                 for j=1:length(UniqueMinIndexParticles)
%                     %Which particle does this spot go to?
%                     CurrentParticleToAssign=ParticleToAssign(UniqueMinIndexParticles(j));
%                     %Which spot needs to be copied onto this particle? I
%                     %need to consider that I might have more than one spot
%                     %that could be assigned to a particle. In this case,
%                     %I'll use distance.
% 
%                     if sum(MinIndexParticles==UniqueMinIndexParticles(j))==1
%                         SpotIndexToCopy=SpotToParticleIndices(find(MinIndexParticles==UniqueMinIndexParticles(j)));
%                     else
%                         [Dummy,SortOrder]=sort(MinValuesParticles(MinIndexParticles==UniqueMinIndexParticles(j)));
%                         SpotIndicesTemp=SpotToParticleIndices(MinIndexParticles==UniqueMinIndexParticles(j));
%                         SpotIndexToCopy=SpotIndicesTemp(SortOrder(1));
%                     end
%                     
%                     %Finally, copy the information onto this particle.
%                     Particles(CurrentParticleToAssign).Index(end+1)=SpotIndexToCopy;
%                     Particles(CurrentParticleToAssign).Frame(end+1)=CurrentFrame;
%                     NewParticlesFlag(SpotIndexToCopy)=0;
%                 end
%             else
%                 %If we have only one previous particle, I need to be more
%                 %careful about the values MinIndexParticles and MinValuesParticles.
%                 
%                 %Which particle does this spot go to?
%                 CurrentParticleToAssign=ParticleToAssign;
%                 %Copy the information onto this particle
%                 Particles(CurrentParticleToAssign).Frame(end+1)=CurrentFrame;
%                 Particles(CurrentParticleToAssign).Index(end+1)=...
%                     SpotToParticleIndices(MinIndexParticles);
%                 NewParticlesFlag(SpotToParticleIndices(MinIndexParticles))=0;
%             end
        end
    end


    %Assign the particles that correspond to a new nucleus
    NewParticlesIndices=find(NewParticlesFlag==1);

    %Indices of the particles that won't be assigned to a nucleus and
    %therefor need to be disapproved by flagging them in SpotFilter
    IndexToMove=[];

   
    for i=1:length(NewParticlesIndices)

        %Recalculate the assigned nuclei in this frame. Note that this is different
        %from AssignedSchnitz in the sense that it asks how many Particles exist in
        %this frame associated with a given Schnitz. 
        %This is useful in case two new particles are closed to a given new nucleus. Right now it
        %will assign the first particle found. I might have to change
        %this if it becomes too annoying.
        AssignedSchnitzCurrentFrame=[];         %These keeps track of which nuclei have
                                    %already been assigned to particles
        for j=1:length(Particles)
            if sum(Particles(j).Frame==CurrentFrame)
                AssignedSchnitzCurrentFrame=[AssignedSchnitzCurrentFrame,Particles(j).Nucleus];
            end
        end


        %Make sure the total number of particles assigned to this schnitz
        %is not higher than SpotsPerNucleus
        if sum(AssignedSchnitzCurrentFrame==MinIndexSchnitz(NewParticlesIndices(i)))<SpotsPerNucleus
            Particles(end+1).Frame=CurrentFrame;
            Particles(end).Index=NewParticlesIndices(i);
            Particles(end).Nucleus=MinIndexSchnitz(NewParticlesIndices(i));
            if retrack
                Particles(end).Approved=0;
            end
        else
            %If the spot cannot be assigned to a particle and schnitz, then we'll
            %take the spot out using SpotFilter.
            IndexToMove=[IndexToMove,i];
        end
    end

    if ~isempty(IndexToMove)     
        SpotFilter(CurrentFrame,IndexToMove)=0;
    end


end