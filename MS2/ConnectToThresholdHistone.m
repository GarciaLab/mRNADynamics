function [fad,fad2,Particles] = ConnectToThresholdHistone(fad,fad2,Particles,...
    CurrentParticle,CurrentFrame,TargetFrame,schnitzcells,Ellipses,SearchRadius,PixelSize)

%V3: Changed to support Laurent's schnitzcells

%This function looks and assigns particle suspects from fad2 into the
%Particles structure.

%First, see which nuclei the fad2 suspects are close to
Particles2=[];

%I'm telling AssignParticle2NucleusV3 to leave any unassigned particles in
%fad2 nevertheless. I hope this is not screwing up the indexing.
[Particles2,fad2,Dummy]=AssignParticle2Nucleus(schnitzcells,Ellipses,Particles2,fad2,fad2,...
    TargetFrame,PixelSize,SearchRadius);


if ~isempty(Particles2)

    %Which one of those corresponds to the nucleus we are looking at?
    NewParticlesNuclei=[Particles2.Nucleus];
    NewParticleSuspect=find(NewParticlesNuclei==Particles(CurrentParticle).Nucleus);

    %Cover the different cases of candidates
    if length(NewParticleSuspect)==1  %If there is only one candidate

        %If there isn't any frame overlap then we can join the traces
        FramesParticle=Particles(CurrentParticle).Frame;
        FramesParticle2=Particles2(NewParticleSuspect).Frame;

        for i=1:length(FramesParticle2)
            FrameOverlap(i,:)=(FramesParticle==FramesParticle2(i));
        end

        if ~sum(sum(FrameOverlap))      %There is no overlap
            SourceFrame=Particles2(NewParticleSuspect).Frame;
            for j=1:length(SourceFrame)
                [fad,fad2,Particles]=...
                    TransferParticle(fad,SourceFrame(j),CurrentParticle,fad2,CurrentFrame,...
                    Particles2(NewParticleSuspect).Index(Particles2(NewParticleSuspect).Frame==SourceFrame(j)),...
                    Particles);
            end
        else
            1+1;    
        end

    elseif isempty(NewParticleSuspect)
        %Do nothing if there are no suspects

    elseif length(NewParticleSuspect)>1
        %If there is more than one particle suspect we need to figure out
        %who's closest to original one
        
        %Find the closest frame in the current particle (the one we'd be
        %copying into)
        
        [MinValue,MinIndex]=min((Particles(CurrentParticle).Frame-TargetFrame).^2);

        %Get the position of the closest current particle (in terms of
        %time)
        [x,y]=fad2xyzFit(Particles(CurrentParticle).Frame(MinIndex),fad, 'addMargin'); 
        xParticle=x(Particles(CurrentParticle).Index(MinIndex));
        yParticle=y(Particles(CurrentParticle).Index(MinIndex));

        %Get the positions of the candidate particles
        Particles2Indices=[Particles2(NewParticleSuspect).Index];
        [x,y]=fad2xyzFit(Particles2(NewParticleSuspect(1)).Frame,fad2, 'addMargin'); 
        xParticles2=x(Particles2Indices);
        yParticles2=y(Particles2Indices);
        
        if length(Particles2(NewParticleSuspect(1)).Index)>1
            1+1;        %I need to think about this case
        end
        
        
        %Now, compare the positions
        clear Distance

        for j=1:length(xParticles2)
            Distance(j,:)=sqrt((xParticles2(j)*PixelSize-...
                xParticle*PixelSize).^2+...
                (yParticles2(j)*PixelSize-...
                yParticle*PixelSize).^2);
        end
        %MinIndex tells us which one of the Particles2 is closest
        [MinValues,MinIndex]=min(Distance');
        
        

        [fad,fad2,Particles]=...
            TransferParticle(fad,Particles2(NewParticleSuspect(1)).Frame,...
            CurrentParticle,fad2,CurrentFrame,...
            Particles2(Particles2Indices(MinIndex)).Index,...
            Particles);
    end

end

% 
% 
% [LastParticleX,LastParticleY]=fad2xyz(CurrentFrame,fad,'addMargin');
% 
% if ~isempty(LastParticleX)
%     if CurrentFrame<TargetFrame
%         LastParticleX=LastParticleX(Particles(CurrentParticle).Index(end));
%         LastParticleY=LastParticleY(Particles(CurrentParticle).Index(end));
%     else
%         LastParticleX=LastParticleX(Particles(CurrentParticle).Index(1));
%         LastParticleY=LastParticleY(Particles(CurrentParticle).Index(1));
%     end
% 
% 
%     [NewParticlesX,NewParticlesY]=fad2xyz(TargetFrame,fad2, 'addMargin');
% 
%     CloseSuspects=sqrt((NewParticlesX-LastParticleX).^2+...
%         (NewParticlesY-LastParticleY).^2)<SearchRadius;
% 
%     if sum(CloseSuspects)
%         if sum(CloseSuspects)>1
%             %If there are many close guys use the one wiht the highest
%             %intensity
%             DogIntensitiesSuspects=fad2.channels(TargetFrame).fits.dog;
%             DogIntensitiesSuspects(~CloseSuspects)=0;
%             [Dummy,MaxSuspectIndex]=max(DogIntensitiesSuspects);
%             CloseSuspects=0;
%             CloseSuspects(MaxSuspectIndex)=1;
%         end
% 
%         [fad,fad2,Particles]=...
%             TransferParticle(fad,TargetFrame,CurrentParticle,fad2,...
%             CurrentFrame,find(CloseSuspects),Particles);
%     end
% end
%     

