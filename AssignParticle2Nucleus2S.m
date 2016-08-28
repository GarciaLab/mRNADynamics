function [Particles,SpotFilter,schnitzcells]=AssignParticle2Nucleus2S(schnitzcells,Ellipses,Particles,Spots,SpotFilter,...
    CurrentFrame,PixelSize,SearchRadius)

%Find which nuclei each of these new particles is closest to.

%V4: Changed to support Laurent's version of schnitzcells and ellipses
%V3: I changed this to be more of an assignment based on the schnitzcells
%tracking.
%V2: I changed this to only detect one particle per nucleus. If there are
%many found it keeps the closest one to the previous particle. Maybe I
%should add some condition related to brightness.


%Assign each particle to a nucleus

%Check if this Particles structure already has the Approved field. If so we
%are performing retracking
if isfield(Particles,'Approved')
    Retracking=1;
    SearchRadius=SearchRadius*2;
else
    Retracking=0;
    SearchRadius=inf;
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
    %MinIndex is a row vector. The position in MinIndex
    %corresponds to each particle. The value at that
    %position corresponds to the closest nucleus.
    
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
    
    [MinValues,MinIndex]=min(Distance');
    
    if Retracking
        MinIndex(MinValues==inf)=inf;
    end


    %Find the schnitz corresponding to MinIndex in the schnitzcell structure.
    %MinIndexNuclei now refers to the identities used in the schnitzcells
    %structure.
    MinIndexNuclei=zeros(size(MinIndex));
    for i=1:length(MinIndex)
        for j=1:length(schnitzcells)
            %I'm not sure this is a good idea, but this is to take care of
            %schnitzses that are inconsistent in their number of cellno
            %entries
            try
                if schnitzcells(j).cellno(schnitzcells(j).frames==(CurrentFrame))==MinIndex(i)
                    MinIndexNuclei(i)=j;
                end
            catch
                display('Nuclear problem')
                %error('What is going on here?')
            end
        end
        if MinIndex(i)==inf
            MinIndexNuclei(i)=inf;
        end
    end

    %Sometimes a nucleus is not in schintzcells. This is probably because
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
    %particle per nucleus. The particles not assigned will moved filtered
    %using SpotFilter.


    AssignedNuclei=[];
    for i=1:length(Particles)
        if ~isempty(Particles(i).Nucleus)
            AssignedNuclei(i)=Particles(i).Nucleus;
        else
            AssignedNuclei(i)=0;
        end
    end


    
    %Find the nuclei in the current frame that were closest to the
    %particles found.
    UniqueMinIndexNuclei=unique(MinIndexNuclei);


    for i=1:length(UniqueMinIndexNuclei)
        
        %Find the particles that are assigned to this nucleus
        ParticlesToAssign=find(AssignedNuclei==UniqueMinIndexNuclei(i));
        %These are the particles that are assigned to the previous nuclei
        
        if ~isempty(ParticlesToAssign)
            for l = 1:length(ParticlesToAssign)
                ParticleToAssign = ParticlesToAssign(l);
                if ~sum(Particles(ParticleToAssign).Frame==CurrentFrame)

                   if (~isempty(ParticleToAssign))&...
                            (sum(MinIndexNuclei==UniqueMinIndexNuclei(i))>=1)
                        %one or more particles are assigned to the same previous
                        %nucleus.


                        %Find the previous particles assigned to this nucleus
                        PreviousParticleIndex=find(AssignedNuclei==UniqueMinIndexNuclei(i));

                        %Get the last positions of the previous particles
                        clear MinDistance
                        clear MinIndex

                        MinDistance = inf(1,length(PreviousParticleIndex));
                        MinIndex = zeros(1,length(PreviousParticleIndex));
                        MinFilter=find(MinIndexNuclei==UniqueMinIndexNuclei(i));

                        for k = 1:length(PreviousParticleIndex)
                            [PreviousParticlesX,PreviousParticlesY]=...
                                SpotsXYZ(Spots(Particles(PreviousParticleIndex(k)).Frame(end)));

                            PreviousParticleX=PreviousParticlesX(Particles(PreviousParticleIndex(k)).Index(end));
                            PreviousParticleY=PreviousParticlesY(Particles(PreviousParticleIndex(k)).Index(end));

                            NewNewParticlesX=NewParticlesX(MinFilter);
                            NewNewParticlesY=NewParticlesY(MinFilter);

                            %Calculate their distances
                            clear Distance

                            for j=1:length(NewNewParticlesX)
                                Distance(j,:)=sqrt((NewNewParticlesX(j)*PixelSize-...
                                    PreviousParticleX*PixelSize).^2+...
                                    (NewNewParticlesY(j)*PixelSize-...
                                    PreviousParticleY*PixelSize).^2);
                            end
                            %Bring the distance of the particles not allowed by
                            %SpotFilter to infinity
                            SpotFilter(CurrentFrame,:);     %HG: Do I need this?

                            %MinIndex tells us which one of the candidates is
                            %closer
                            [MinDistance(k),MinIndex(k)]=min(Distance');
                        end

                        %Case where two particles are different
                        if length(unique(MinFilter(MinIndex)))==length(MinFilter(MinIndex))
                            for  k = 1:length(PreviousParticleIndex)
                            %This is the index of the particle found
                                if ((MinDistance(k))<SearchRadius*PixelSize)&...
                                        ~sum(Particles(PreviousParticleIndex(k)).Frame(end)==CurrentFrame)
                                    Particles(PreviousParticleIndex(k)).Frame(end+1)=CurrentFrame;
                                    Particles(PreviousParticleIndex(k)).Index(end+1)=MinFilter(MinIndex(k));
                                    NewParticlesFlag(MinFilter(MinIndex(k)))=0;
                                end
                            end

                        %case where they detected the same particle
                        else
                            [CloserDistance,CloserIndex]=min(MinDistance);
                            %This is the index of the particle found
                                if ((MinDistance(CloserIndex))<SearchRadius*PixelSize)&...
                                        ~sum(Particles(PreviousParticleIndex(CloserIndex)).Frame(end)==CurrentFrame)
                                    Particles(PreviousParticleIndex(CloserIndex)).Frame(end+1)=CurrentFrame;
                                    Particles(PreviousParticleIndex(CloserIndex)).Index(end+1)=MinFilter(MinIndex(CloserIndex));
                                    NewParticlesFlag(MinFilter(MinIndex(CloserIndex)))=0;
                                end  
                        end
                    end
                end
            end
        end
    end


    %Assign the particles that correspond to a new nucleus
    NewParticlesIndices=find(NewParticlesFlag==1);

    %Indices of the particles that need to be filtered out using
    %SpotFilter.
    IndexToMove=[];

    for i=1:length(NewParticlesIndices)

        %Recalculate the assigned nuclei. This is useful in case two
        %new particles are close to a given new nucleus. Right now it
        %will assign the first particle found. I might have to change
        %this if it becomes too annoying.
        AssignedNuclei=[];          %These keeps track of which nuclei have
                                    %already been assigned to particles
        AssignedParticles=[];
        for j=1:length(Particles)
            AssignedNuclei=[AssignedNuclei,Particles(j).Nucleus];
            %AssignedParticles=[AssignedParticles,Particles(j).Index(end)];
        end


        %Make sure this new particle doesn't get assigned to an
        %already-assigned nucleus (Need to modify this to be okay for 2
        %spots per nucleus)
        if (sum(AssignedNuclei==MinIndexNuclei(NewParticlesIndices(i)))<=1)
            Particles(end+1).Frame=CurrentFrame;
            Particles(end).Index=NewParticlesIndices(i);
            Particles(end).Nucleus=MinIndexNuclei(NewParticlesIndices(i));
            if Retracking
                Particles(end).Approved=0;
            end
            
        %If the particle can't be put anywhere else then move it to
        %fad2. In order to do this we will create a list of the indices
        %that need to be moved. This is important so that we don't
        %lose track of anything.
        else
            IndexToMove=[IndexToMove,i];
        end
    end

    %Now use the SpotFilter array to stop considering the particles that were not assigned
    if ~isempty(IndexToMove)
        
        SpotFilter(CurrentFrame,IndexToMove)=0;

        %Since we got rid of particles in fad we need to shift the
        %index entry for the remaining particles accordingly
        for i=1:length(Particles)
            if sum(Particles(i).Frame==CurrentFrame)
                Particles(i).Index(Particles(i).Frame==CurrentFrame)=...
                    Particles(i).Index(Particles(i).Frame==CurrentFrame)-...
                    sum(NewParticlesIndices(IndexToMove)<Particles(i).Index(Particles(i).Frame==CurrentFrame));
            end
        end





    end


end