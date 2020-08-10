function ParticleOutput = identifyParticle(Spots, Particles, CurrentFrame, ...
    CurrentChannelIndex, UseHistoneOverlay, schnitzcells, ConnectPosition)
%IDENTIFYPARTICLE Summary of this function goes here
%   Detailed explanation goes here

numParticles = length(Particles{CurrentChannelIndex});
if ~exist('ConnectPosition', 'var')
    [ConnectPositionx,ConnectPositiony]=ginput(1);
    ConnectPosition = [ConnectPositionx,ConnectPositiony];
end

if ~isempty(ConnectPosition)
    
    display(ConnectPosition);
    
    %Find the closest particle
    [ParticleOutput,~]=FindClickedParticle(ConnectPosition,CurrentFrame,...
        Spots{CurrentChannelIndex},Particles{CurrentChannelIndex});
    
    disp(['Clicked particle: ',num2str(ParticleOutput)]);

    if UseHistoneOverlay
        %Find the closest nucleus
        NewNuclei=ConnectPosition;

        %Find which schnitz this corresponds to
        SchnitzSuspect=[];
        xPosSuspect=[];
        yPosSuspect=[];
        for j=1:length(schnitzcells)
            if sum(schnitzcells(j).frames==CurrentFrame)
                SchnitzSuspect=[SchnitzSuspect,j];
                xPosSuspect=[xPosSuspect,...
                    schnitzcells(j).cenx(find((schnitzcells(j).frames)==CurrentFrame))];
                yPosSuspect=[yPosSuspect,...
                    schnitzcells(j).ceny(find((schnitzcells(j).frames)==CurrentFrame))];
            end
        end

        %Find the closest one to the point where we clicked
        Distance=sqrt(double((NewNuclei(1)-xPosSuspect).^2+(NewNuclei(2)-yPosSuspect).^2));
        [~,ClosestNucleusIndex]=min(Distance);

        ClickedSchnitz=SchnitzSuspect(ClosestNucleusIndex);

        %Now, find its associated particle
        for i=1:numParticles
            if ~isempty(Particles{CurrentChannelIndex}(i).Nucleus)
                AssignedNuclei(i)=Particles{CurrentChannelIndex}(i).Nucleus;
            else
                AssignedNuclei(i)=nan;
            end
        end
        AssociatedParticle=find(AssignedNuclei==ClickedSchnitz);

        if isempty(AssociatedParticle)
            disp(['Nucleus ',num2str(ClickedSchnitz),' does not have an associated particle'])
        else
            disp(['Particle ',num2str(AssociatedParticle),' is associate with nucleus ',...
                num2str(ClickedSchnitz)])
        end
    end
end
end

