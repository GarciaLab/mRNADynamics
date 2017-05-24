function [Particles,SpotFilter,schnitzcells]=AssignParticle2NucleusV2(...
    schnitzcells,Ellipses,Particles,Spots,SpotFilter,...
    CurrentFrame,PixelSize,SearchRadius,SpotsPerNucleus)






%Find the N closest particles to each nucleus. N is determined by
%ExperimentType (1spot, 2spots). If N>1, then the code also uses the spot
%position for tracking.

%SearchRadius is given in pixels, but it comes from a setting in um
%from TrackmRNADynamics.

%Check if this Particles structure already has the Approved field. If so we
%are performing retracking
if isfield(Particles,'Approved')
    Retracking=1;
    %HG for v2, I got rid of this
    %SearchRadius=SearchRadius*2;        %HG: Why did I decide to do this?
else
    Retracking=0;
    %SearchRadius=inf;
end


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
    clear Distance

    for j=1:length(NewSpotsX)
        Distance(j,:)=sqrt((NewSpotsX(j)*PixelSize-...
            NewEllipsesX*PixelSize).^2+...
            (NewSpotsY(j)*PixelSize-...
            NewEllipsesY*PixelSize).^2);
    end
        
    %If retracking make the distance for the particles that have been
    %approved infinite.
    if Retracking
        for i=1:length(Particles)
            %Find which approved particles are in this frame
            if sum(Particles(i).Frame==CurrentFrame)&(Particles(i).Approved~=0)
                %Make the distance infinite
                Distance(Particles(i).Index(find(Particles(i).Frame==CurrentFrame)),:)=inf;
                NewParticlesFlag(Particles(i).Index(find(Particles(i).Frame==CurrentFrame)))=-1;
            end
        end
    end
    
%     %Set any distances that are larger than SearchRadius to infinity
%     Distance(Distance>SearchRadius)=inf;    
    
    %MinIndex is a row vector. The position in MinIndex
    %corresponds to each spot. The value at that
    %position corresponds to the closest ellipse.
    [MinValues,MinIndex]=min(Distance');
    
%     %Allow only the closest SpotsPerNucleus spots to be close to each
%     %Ellipse.
%     %Find the unique Ellipses
%     UniqueEllipses=unique(MinIndex);
%     %Check how many spots are associated with each UniqueEllipse. If there
%     %are more than SpotsPerNucleus spots closest to a given Ellipse, we'll
%     %set their distances to infinity and recalculate MinValues and
%     %MinIndex.
%     
%     %How many spots are closest to a nucleus?
%     NSpotsPerNucleusCurrent=zeros(size(UniqueEllipses));
%     for j=1:length(UniqueEllipses)
%         NSpotsPerNucleusCurrent(j)=sum(MinIndex==UniqueEllipses(j));
%     end
%     
%     MinValues(MinIndex==17)
%     
%     for j=1:length(UniqueEllipses)
%         if sum(MinIndex==UniqueEllipses(j))>SpotsPerNucleus
%             %Find the indices of teh Spots that are closest to this Ellipse
%             %from MinIndex.
%             ClosestSpots=find(MinIndex==UniqueEllipses(j));
%             %Sort the Spots according to their distance to the Ellipse.
%             [SortedValues,SortedIndices]=sort(MinValues(MinIndex==UniqueEllipses(j)));
%             %Set the distance of the remaining Spots to this Ellipse to
%             %infinity
%             Distance(ClosestSpots(SortedIndices(SpotsPerNucleus+1:end)),UniqueEllipses(j))=inf;
%         end
%     end
%     
%     %Recalculate MinValues and MinIndex in case we made some changes in the
%     %for-loop above.
%     [MinValues,MinIndex]=min(Distance');
%     %Recalculate the unique ellipses
%     UniqueEllipses=unique(MinIndex);
%     %Recalculate the number of spots closest to each nucleus
%     NSpotsPerNucleusCurrent=zeros(size(UniqueEllipses));
%     for j=1:length(UniqueEllipses)
%         NSpotsPerNucleusCurrent(j)=sum(MinIndex==UniqueEllipses(j));
%     end
%     
%     
%     
%     while sum(NSpotsPerNucleusCurrent>SpotsPerNucleus)
%         for j=1:length(UniqueEllipses)
%             if sum(MinIndex==UniqueEllipses(j))>SpotsPerNucleus
%                 %Find the indices of teh Spots that are closest to this Ellipse
%                 %from MinIndex.
%                 ClosestSpots=find(MinIndex==UniqueEllipses(j));
%                 %Sort the Spots according to their distance to the Ellipse.
%                 [SortedValues,SortedIndices]=sort(MinValues(MinIndex==UniqueEllipses(j)));
%                 %Set the distance of the remaining Spots to this Ellipse to
%                 %infinity
%                 Distance(ClosestSpots(SortedIndices(SpotsPerNucleus+1:end)),UniqueEllipses(j))=inf;
%             end
%         end
%         %Recalculate MinValues and MinIndex in case we made some changes in the
%         %for-loop above.
%         [MinValues,MinIndex]=min(Distance');
%         %Recalculate the unique ellipses
%         UniqueEllipses=unique(MinIndex);
%         %Recalculate the number of spots closest to each nucleus
%         NSpotsPerNucleusCurrent=zeros(size(UniqueEllipses));
%         for j=1:length(UniqueEllipses)
%             NSpotsPerNucleusCurrent(j)=sum(MinIndex==UniqueEllipses(j));
%         end
%     end
%         
        
               
    %If retracking, set the distance of the already assigned spots to
    %infinity
    if Retracking
        MinIndex(MinValues==inf)=inf;
    end
    
    
    %Find the schnitz corresponding to MinIndex in the schnitzcell structure.
    %MinIndexNuclei is now a row vector. The position in MinIndexNuclei
    %corresponds to each spot. The value at that position corresponds
    %to the closest schnitz.
    MinIndexSchnitz=zeros(size(MinIndex));
    for i=1:length(MinIndex)
        for j=1:length(schnitzcells)
            if schnitzcells(j).cellno(schnitzcells(j).frames==(CurrentFrame))==MinIndex(i)
                MinIndexSchnitz(i)=j;
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

    
    %For debugging
    if CurrentFrame==2
        1+1;
    end
    
    
    
    for i=1:length(UniqueMinIndexSchnitz)
        
        %Find the particles that are already assigned to this schnitz
        ParticleToAssign=find(AssignedSchnitz==UniqueMinIndexSchnitz(i));

        %Are there any particles in previous frames that are assigned to
        %this schnitz? If not, we move on and define the current spots as
        %new particles.
        if ~isempty(ParticleToAssign)
            %Spots I need to locate to this schnitz
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


            %Calculate the distances between the spots on this frame and the
            %ellipses
            clear DistanceParticles

            for j=1:length(SpotToParticleX)
                DistanceParticles(j,:)=sqrt((SpotToParticleX(j)*PixelSize-...
                    PreviousParticlesX*PixelSize).^2+...
                    (SpotToParticleY(j)*PixelSize-...
                    PreviousParticlesY*PixelSize).^2);
            end
            %MinIndexParticles is a row vector. The position in MinIndex
            %corresponds to each spot. The value at that
            %position corresponds to the closest particle. Note that this
            %is not true if we have only one particle
            [MinValuesParticles,MinIndexParticles]=min(DistanceParticles');
            
            
            %Sometimes the two or more new Spots are closer to the same
            %previous Particle. In that case, only keep the closest
            %Spot-Particle pair and set the distance of the other spots to
            %the same particle to infinity.
            UniqueMinIndexParticles=unique(MinIndexParticles);
            
                        
            %Assign the spots to their corresponding particles. I need to
            %consider the case where there was only one previous particle.
            if length(ParticleToAssign)>1
                %If we have more than one previous particle, I can make use
                %of the matrix Distance, of MinValuesParticles, and
                %MinIndexParticles.
                for j=1:length(UniqueMinIndexParticles)
                    %Which particle does this spot go to?
                    CurrentParticleToAssign=ParticleToAssign(UniqueMinIndexParticles(j));
                    %Which spot needs to be copied onto this particle? I
                    %need to consider that I might have more than one spot
                    %that could be assigned to a particle. In this case,
                    %I'll use distance.

                    if sum(MinIndexParticles==UniqueMinIndexParticles(j))==1
                        SpotIndexToCopy=SpotToParticleIndices(find(MinIndexParticles==UniqueMinIndexParticles(j)));
                    else
                        [Dummy,SortOrder]=sort(MinValuesParticles(MinIndexParticles==UniqueMinIndexParticles(j)));
                        SpotIndicesTemp=SpotToParticleIndices(MinIndexParticles==UniqueMinIndexParticles(j));
                        SpotIndexToCopy=SpotIndicesTemp(SortOrder(1));
                    end
                    
                    %Finally, copy the information onto this particle.
                    Particles(CurrentParticleToAssign).Index(end+1)=SpotIndexToCopy;
                    Particles(CurrentParticleToAssign).Frame(end+1)=CurrentFrame;
                    NewParticlesFlag(SpotIndexToCopy)=0;
                end
            else
                %If we have only one previous particle, I need to be more
                %careful about the values MinIndexParticles and MinValuesParticles.
                
                %Which particle does this spot go to?
                CurrentParticleToAssign=ParticleToAssign;
                %Copy the information onto this particle
                Particles(CurrentParticleToAssign).Frame(end+1)=CurrentFrame;
                Particles(CurrentParticleToAssign).Index(end+1)=...
                    SpotToParticleIndices(MinIndexParticles);
                NewParticlesFlag(SpotToParticleIndices(MinIndexParticles))=0;
            end
            
%             
%             
%             for j=1:length(UniqueMinIndexParticles)
%                 if sum(MinIndexParticles==UniqueMinIndexParticles(j))>1
%                     %Find the indices of teh Spots that are closest to this Ellipse
%                     %from MinIndex.
%                     ClosestSpots=find(MinIndexParticles==UniqueMinIndexParticles(j));
%                     %Sort the Spots according to their distance to the Ellipse.
%                     [SortedValues,SortedIndices]=...
%                         sort(MinValuesParticles(MinIndexParticles==UniqueMinIndexParticles(j)));
%                     %Set the distance of the remaining Spots to this Ellipse to
%                     %infinity
%                     DistanceParticles(ClosestSpots(SortedIndices(2:end)),UniqueMinIndexParticles(j))=inf;
%                 end
%             end

%             
%             
%             
%             %I need to keep doing this until            
%             if length(UniqueMinIndexParticles)<length(MinIndexParticles)
%                 for j=1:length(UniqueMinIndexParticles)
%                     %Check whether more than one spot is matched to the same
%                     %particle.
%                     if sum(MinIndexParticles==UniqueMinIndexParticles(j))>1
%                         [Dummy,SortOrder]=...
%                             sort(MinValuesParticles(MinIndexParticles==UniqueMinIndexParticles(j)));
%                         DistanceParticles(SortOrder(2:end),UniqueMinIndexParticles(j))=inf;
%                     end
%                 end
%                 [MinValuesParticles,MinIndexParticles]=min(DistanceParticles');
%                 UniqueMinIndexParticles=unique(MinIndexParticles);
%             end
            
%             
%             %I need to keep doing this until            
%             while length(UniqueMinIndexParticles)<length(MinIndexParticles)
%                 for j=1:length(UniqueMinIndexParticles)
%                     %Check whether more than one spot is matched to the same
%                     %particle.
%                     if sum(MinIndexParticles==UniqueMinIndexParticles(j))>1
%                         [Dummy,SortOrder]=...
%                             sort(MinValuesParticles(MinIndexParticles==UniqueMinIndexParticles(j)));
%                         DistanceParticles(SortOrder(2:end),UniqueMinIndexParticles(j))=inf;
%                     end
%                 end
%                 [MinValuesParticles,MinIndexParticles]=min(DistanceParticles');
%                 UniqueMinIndexParticles=unique(MinIndexParticles);
%             end
%             
%             
% 
%             %Assign the spots to their corresponding particles. I need to
%             %consider the case where there was only one previous particle.
%             if length(ParticleToAssign)>1
%                 %If we have more than one previous particle, I can make use
%                 %of the matrix Distance, of MinValuesParticles, and
%                 %MinIndexParticles.
%                 for j=1:length(MinIndexParticles)
%                     %Which particle does this spot go to?
%                     CurrentParticleToAssign=ParticleToAssign(MinIndexParticles(j));
%                     %Copy the information onto this particle
%                     Particles(CurrentParticleToAssign).Frame(end+1)=CurrentFrame;
%                     Particles(CurrentParticleToAssign).Index(end+1)=SpotToParticleIndices(j);
%                     NewParticlesFlag(SpotToParticleIndices(j))=0;
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
% 
% 
%             
            
            
                        
%             if ~sum(Particles(ParticleToAssign).Frame==CurrentFrame)
%                 if ((sum(MinIndexSchnitz==UniqueMinIndexSchnitz(i)))==1)&(~isempty(ParticleToAssign))&...
%                         (MinValues(MinIndexSchnitz==UniqueMinIndexSchnitz(i)))<SearchRadius*PixelSize
%                     %One particle is assigned to a previous one within a nucleus
% 
%                     Particles(ParticleToAssign).Frame(end+1)=CurrentFrame;
%                     Particles(ParticleToAssign).Index(end+1)=find(MinIndexSchnitz==...
%                                     UniqueMinIndexSchnitz(i));
%                     NewParticlesFlag(find(MinIndexSchnitz==...
%                         UniqueMinIndexSchnitz(i)))=0;
% 
% 
%                 elseif (~isempty(ParticleToAssign))&...
%                         (sum(MinIndexSchnitz==UniqueMinIndexSchnitz(i))>1)
%                     %Two or more particles are assigned to the same previous
%                     %nucleus.
% 
% 
%                     %Find the previous particle assigned to this nucleus
%                     PreviousParticleIndex=find(AssignedNuclei==UniqueMinIndexNuclei(i));
% 
%                     %Get the last position of the previous particle
%                     [PreviousParticlesX,PreviousParticlesY]=...
%                         SpotsXYZ(Spots(Particles(PreviousParticleIndex).Frame(end)));
%                                                         
%                     PreviousParticleX=PreviousParticlesX(Particles(PreviousParticleIndex).Index(end));
%                     PreviousParticleY=PreviousParticlesY(Particles(PreviousParticleIndex).Index(end));
% 
%                     MinFilter=find(MinIndexNuclei==UniqueMinIndexNuclei(i));
% 
%                     NewNewSpotsX=NewSpotsX(MinFilter);
%                     NewNewSpotsY=NewSpotsY(MinFilter);
% 
%                     %Calculate their distances
%                     clear Distance
% 
%                     for j=1:length(NewNewSpotsX)
%                         Distance(j,:)=sqrt((NewNewSpotsX(j)*PixelSize-...
%                             PreviousParticleX*PixelSize).^2+...
%                             (NewNewSpotsY(j)*PixelSize-...
%                             PreviousParticleY*PixelSize).^2);
%                     end
%                     %MinIndex tells us which one of the candidates is
%                     %closer
%                     [MinValues2,MinIndex]=min(Distance');
% 
%                     %This is the index of the particle found
%                     if Distance(MinIndex)<SearchRadius*PixelSize
%                         Particles(PreviousParticleIndex).Frame(end+1)=CurrentFrame;
%                         Particles(PreviousParticleIndex).Index(end+1)=MinFilter(MinIndex);
%                         NewParticlesFlag(MinFilter(MinIndex))=0;
%                     end
%                 end
%             end
        end
    end


    %Assign the particles that correspond to a new nucleus
    NewParticlesIndices=find(NewParticlesFlag==1);

    %Indices of the particles that won't be assigned to a nucleus and
    %therefor need to be disapproved by flagging them in SpotFilter
    IndexToMove=[];

    for i=1:length(NewParticlesIndices)

        %Recalculate the assigned nuclei. This is useful in case two
        %new particles are closed to a given new nucleus. Right now it
        %will assign the first particle found. I might have to change
        %this if it becomes too annoying.
        AssignedSchnitz=[];         %These keeps track of which nuclei have
                                    %already been assigned to particles
        for j=1:length(Particles)
            AssignedSchnitz=[AssignedSchnitz,Particles(j).Nucleus];
        end


        %Make sure the total number of particles assigned to this schnitz
        %is not higher than SpotsPerNucleus
        if sum(AssignedSchnitz==MinIndexSchnitz(NewParticlesIndices(i)))<SpotsPerNucleus
            Particles(end+1).Frame=CurrentFrame;
            Particles(end).Index=NewParticlesIndices(i);
            Particles(end).Nucleus=MinIndexSchnitz(NewParticlesIndices(i));
            if Retracking
                Particles(end).Approved=0;
            end

        %If the spot cannot be assigned to a particle and schnitz, then we'll
        %take the spot out using SpotFilter.
            IndexToMove=[IndexToMove,i];
        end
    end

    if ~isempty(IndexToMove)     
        SpotFilter(CurrentFrame,IndexToMove)=0;
    end


end