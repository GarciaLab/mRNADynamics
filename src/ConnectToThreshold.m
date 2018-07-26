function [fad,fad2,Particles] = ConnectToThreshold2(fad,fad2,Particles,...
    CurrentParticle,CurrentFrame,TargetFrame,SearchRadius)

%This function looks and assigns particle suspects from fad2 into the
%Particles structure.

[LastParticleX,LastParticleY]=fad2xyz(CurrentFrame,fad,'addMargin');

if ~isempty(LastParticleX)
    if CurrentFrame<TargetFrame
        LastParticleX=LastParticleX(Particles(CurrentParticle).Index(end));
        LastParticleY=LastParticleY(Particles(CurrentParticle).Index(end));
    else
        LastParticleX=LastParticleX(Particles(CurrentParticle).Index(1));
        LastParticleY=LastParticleY(Particles(CurrentParticle).Index(1));
    end


    [NewParticlesX,NewParticlesY]=fad2xyz(TargetFrame,fad2, 'addMargin');

    CloseSuspects=sqrt((NewParticlesX-LastParticleX).^2+...
        (NewParticlesY-LastParticleY).^2)<SearchRadius;

    if sum(CloseSuspects)
        if sum(CloseSuspects)>1
            %If there are many close guys use the one with the highest
            %intensity
            DogIntensitiesSuspects=fad2.channels(TargetFrame).fits.dog;
            DogIntensitiesSuspects(~CloseSuspects)=0;
            [Dummy,MaxSuspectIndex]=max(DogIntensitiesSuspects);
            CloseSuspects=0;
            CloseSuspects(MaxSuspectIndex)=1;
        end

        [fad,fad2,Particles]=...
            TransferParticle(fad,TargetFrame,CurrentParticle,fad2,...
            CurrentFrame,find(CloseSuspects),Particles);
    end
end
    

