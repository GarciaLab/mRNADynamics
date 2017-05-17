function [Particles,SpotFilter,schnitzcells]=AssignParticle2Nucleus2S(...
    schnitzcells,Ellipses,Particles,Spots,SpotFilter,...
    CurrentFrame,PixelSize,SearchRadius)

%Find which nuclei each of these new particles is closest to. This code
%supports up to two particles per nucleus.

%Assign each particle to a nucleus

%Check if this Particles structure already has the Approved field. If so we
%are performing retracking
if isfield(Particles,'Approved')
    Retracking=1;
    SearchRadius=SearchRadius*2;
else
    Retracking=0;
    SearchRadius=inf;
    SearchRadiusTwoSpot=2;      %In um. This is the distance used to track
                                %two spots within the same nucleus.
end


%Get the particles and nuclei positions. The order of both the nuclei and
%particle position is raw, not the one in the Particles and Nuclei
%structures.
[NewParticlesX,NewParticlesY]=SpotsXYZ(Spots(CurrentFrame));

if ~isempty(NewParticlesX)

    %Get the position of the nuclei in this frame
    NewNucleiXY=Ellipses{CurrentFrame}(:,[1,2]);
    NewNucleiX=NewNucleiXY(:,1);
    NewNucleiY=NewNucleiXY(:,2);

    %This is to keep track of already-assigned particles
    NewParticlesFlag=ones(size(NewParticlesX));


    %Calculate the distances between the particles on this frame and the
    %nuclei
    clear Distance

    for j=1:length(NewParticlesX)
        Distance(j,:)=sqrt((NewParticlesX(j)*PixelSize-...
            NewNucleiX*PixelSize).^2+...
            (NewParticlesY(j)*PixelSize-...
            NewNucleiY*PixelSize).^2);
    end

    %If a particle already exists in this frame, then take it out of the
    %pool for consideration.
    for i=1:length(Particles)
        %Find which approved particles are in this frame
        if sum(Particles(i).Frame==CurrentFrame)
            %Make the distance infinite
            Distance(Particles(i).Index(find(Particles(i).Frame==CurrentFrame)),:)=inf;
        end
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
    
    %MinIndex is a row vector. The position in MinIndex
    %corresponds to each particle. The value at that
    %position corresponds to the closest ellipse.
    [MinValues,MinIndex]=min(Distance');
    
    if Retracking
        MinIndex(MinValues==inf)=inf;
    end


    %Find the schnitz corresponding to the ellipse given by MinIndex.
    %MinIndexNuclei now refers to the identities used in the schnitzcells
    %structure.
    MinIndexNuclei=zeros(size(MinIndex));
    for i=1:length(MinIndex)
        for j=1:length(schnitzcells)
            if schnitzcells(j).cellno(schnitzcells(j).frames==(CurrentFrame))==MinIndex(i)
                MinIndexNuclei(i)=j;
            end
        end
        if MinIndex(i)==inf
            MinIndexNuclei(i)=inf;
        end
    end

    %Sometimes a nucleus is not in schnitzcells. This is probably because
    %of the user having deleted it manually. In that case, re-create the
    %schnitz.
    
    if sum(find(MinIndexNuclei==0))
        error('Check this. Did schnitzcells get edited?')
        NotAssignedNuclei=find(MinIndexNuclei==0);
        for i=1:length(NotAssignedNuclei)
            i;
            schnitzcells(end+1).P=0;
            schnitzcells(end).E=0;
            schnitzcells(end).D=0;
            schnitzcells(end).frames=CurrentFrame;
            schnitzcells(end).cenx=Ellipses{CurrentFrame}(MinIndex(NotAssignedNuclei(i)),1);
            schnitzcells(end).ceny=Ellipses{CurrentFrame}(MinIndex(NotAssignedNuclei(i)),2);
            schnitzcells(end).ceny=Ellipses{CurrentFrame}(MinIndex(NotAssignedNuclei(i)),5);        
            schnitzcells(end).len=max([Ellipses{CurrentFrame}(MinIndex(NotAssignedNuclei(i)),3),...
                Ellipses{CurrentFrame}(MinIndex(NotAssignedNuclei(i)),4)]);  
            schnitzcells(end).cellno=MinIndex(NotAssignedNuclei(i));
            %schnitzcells(end).approved=0;
            MinIndexNuclei(NotAssignedNuclei(i))=length(schnitzcells);
        end
    end
        


    
    
    %We need to go through all possible cases:
    %1) One particle is assigned to a previous one within a nucleus
    %2) If there are multiple particles closest to a nucleus we pick
    %the closest/brightest one. As a result, there can only be on
    %particle per nucleus. The particles not assigned will disapproved
    %by flagging them in SpotsFilter.
    AssignedNuclei=[];
    for i=1:length(Particles)
        if ~isempty(Particles(i).Nucleus)
            AssignedNuclei(i)=Particles(i).Nucleus;
        else
            AssignedNuclei(i)=0;
        end
    end

    
    %Find the schnitz in the current frame that were closest to the
    %particles found.
    UniqueMinIndexNuclei=unique(MinIndexNuclei);

    for i=1:length(UniqueMinIndexNuclei)
        
        %Find the particle that is assigned to this nucleus
        ParticleToAssign=find(AssignedNuclei==UniqueMinIndexNuclei(i));

        if ~isempty(ParticleToAssign)
        %HG: I removed the following if-statement because I now check above
        %whether particles were arelady assigned to this frame and modify
        %the Distance matrix accordingly
        %if ~sum(Particles(ParticleToAssign).Frame==CurrentFrame)
            if ((sum(MinIndexNuclei==UniqueMinIndexNuclei(i)))==1)&...
                    (MinValues(MinIndexNuclei==UniqueMinIndexNuclei(i)))<SearchRadius*PixelSize
                %One particle is assigned to a previous one within a
                %nucleus.

                Particles(ParticleToAssign).Frame(end+1)=CurrentFrame;
                Particles(ParticleToAssign).Index(end+1)=find(MinIndexNuclei==...
                                UniqueMinIndexNuclei(i));
                NewParticlesFlag(find(MinIndexNuclei==...
                    UniqueMinIndexNuclei(i)))=0;


            elseif sum(MinIndexNuclei==UniqueMinIndexNuclei(i))==2
                %Two particles are assigned to the same previous
                %nucleus. We need to find which new spot
                %is closer to which previous particle.

                %Find the previous particles assigned to this nucleus
                PreviousParticleIndex=find(AssignedNuclei==UniqueMinIndexNuclei(i));

                if length(PreviousParticleIndex)==2
                    1+1
                end

                %Calculate the distances between the spots on this frame and the
                %particles in the previous frame
                
                %Get the position of all particles in the previous frame
                [PreviousParticlesX,PreviousParticlesY]=...
                    SpotsXYZ(Spots(Particles(PreviousParticleIndex(end)).Frame(end)));
                %Filter for the particles in the previous frame that are
                %closer to the new spots in the current frame
                PreviousParticleX=[];
                PreviousParticleY=[];
                for j=1:length(PreviousParticleIndex)
                    PreviousParticleX(j)=PreviousParticlesX(Particles(PreviousParticleIndex(j)).Index(end));
                    PreviousParticleY(j)=PreviousParticlesY(Particles(PreviousParticleIndex(j)).Index(end));
                end
                %Get the two new spots we're going to compare to
                %the previous particles
                MinFilter=find(MinIndexNuclei==UniqueMinIndexNuclei(i));

                NewNewParticlesX=NewParticlesX(MinFilter);
                NewNewParticlesY=NewParticlesY(MinFilter);

                clear Distance              
                for j=1:length(NewNewParticlesX)
                    Distance(j,:)=sqrt((NewNewParticlesX(j)*PixelSize-...
                        PreviousParticleX'*PixelSize).^2+...
                        (NewNewParticlesY(j)*PixelSize-...
                        PreviousParticleY'*PixelSize).^2);
                end
                %MinIndex2Spots is a row vector. The position in MinIndex2Spots
                %corresponds to each new spot. The value at that
                %position corresponds to the closest previous particle.
                [MinValues2Spots,MinIndex2Spots]=min(Distance');

                if length(unique(MinIndex2Spots))~=length(MinIndex2Spots)
                   error('HG should check this case') 
                end
                
                %Go through each distance, make sure that it's within the
                %search radius, and assign the new spots to the previous
                %particles.
                for j=1:length(MinIndex2Spots)
                    if MinValues2Spots(j)<SearchRadiusTwoSpot
                        Particles(PreviousParticleIndex(MinIndex2Spots(j))).Frame(end+1)=CurrentFrame;
                        Particles(PreviousParticleIndex(MinIndex2Spots(j))).Index(end+1)=MinFilter(j);
                        NewParticlesFlag(MinFilter(j))=0;
                    end
                end
%                 
%                 
%                 
%                 for j=1:length(PreviousParticleIndex)
%                     %Get the last position of the previous particle
%                     [PreviousParticlesX,PreviousParticlesY]=...
%                         SpotsXYZ(Spots(Particles(PreviousParticleIndex(j)).Frame(end)));
%                     PreviousParticleX=PreviousParticlesX(Particles(PreviousParticleIndex(j)).Index(end));
%                     PreviousParticleY=PreviousParticlesY(Particles(PreviousParticleIndex(j)).Index(end));
% 
%                     %Get the two new spots we're going to compare to
%                     %the previous particles
%                     MinFilter=find(MinIndexNuclei==UniqueMinIndexNuclei(i));
%                     MinFilter=MinFilter;
% 
%                     NewNewParticlesX=NewParticlesX(MinFilter);
%                     NewNewParticlesY=NewParticlesY(MinFilter);
% 
%                     %Calculate their distances
%                     clear Distance
% 
%                     for k=1:length(NewNewParticlesX)
%                         Distance(k,:)=sqrt((NewNewParticlesX(k)*PixelSize-...
%                             PreviousParticleX*PixelSize).^2+...
%                             (NewNewParticlesY(k)*PixelSize-...
%                             PreviousParticleY*PixelSize).^2);
%                     end
%                     %If the minimum particle has already been assigned,
%                     %then I need to set its distance to infinity
%                     Distance(find(NewParticlesFlag(MinFilter)==0))=inf;
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

            elseif (~isempty(ParticleToAssign))&...
                (sum(MinIndexNuclei==UniqueMinIndexNuclei(i))>2)

                warning('HG still needs to check this')

            end
            %end
        end
    end


    %Assign the particles that correspond to a new nucleus
    NewParticlesIndices=find(NewParticlesFlag==1);

    %Indices of the particles that won't be assigned to a nucleus and
    %therefor need to be disapproved by flagging them in SpotFilter
    IndexToMove=[];

    for i=1:length(NewParticlesIndices)

        %Recalculate the assigned nuclei. This is useful in case two
        %new particles are closed to a given new nucleus. 
        AssignedNuclei=[];          %These keeps track of which nuclei have
                                    %already been assigned to particles
        for j=1:length(Particles)
            AssignedNuclei=[AssignedNuclei,Particles(j).Nucleus];
        end


        %Is this spot close to an already-assigned nucleus?
        
        %If the spot is close to an unassigned nucleus, or if the nucleus
        %already has one other assigned particle, then create a new
        %particle.
        if sum(AssignedNuclei==MinIndexNuclei(NewParticlesIndices(i)))<2
            Particles(end+1).Frame=CurrentFrame;
            Particles(end).Index=NewParticlesIndices(i);
            Particles(end).Nucleus=MinIndexNuclei(NewParticlesIndices(i));
            if Retracking
                Particles(end).Approved=0;
            end
            
        %If the particle can't be put anywhere else then move it out of the
        %pool by taking it out of SpotFilter. This is important so that we 
        %can bring it back if necessary
        else
            IndexToMove=[IndexToMove,i];
        end
    end

    if ~isempty(IndexToMove)     
        SpotFilter(CurrentFrame,IndexToMove)=0;
    end


end