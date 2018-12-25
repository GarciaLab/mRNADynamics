function ParticleOutput = identifyParticle(Spots, Particles, CurrentFrame, ...
    CurrentChannel, UseHistoneOverlay, schnitzcells)
%IDENTIFYPARTICLE Summary of this function goes here
%   Detailed explanation goes here

numParticles = length(Particles{CurrentChannel});
[ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'b', 'linewidth',1);
ConnectPosition = [ConnectPositionx,ConnectPositiony];
if ~isempty(ConnectPosition)
    %Find the closest particle
    [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,Spots{CurrentChannel},Particles{CurrentChannel});
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
        Distance=sqrt((NewNuclei(1)-xPosSuspect).^2+(NewNuclei(2)-yPosSuspect).^2);
        [MinValue,ClosestNucleusIndex]=min(Distance);

        ClickedSchnitz=SchnitzSuspect(ClosestNucleusIndex);

        %Now, find its associated particle
        for i=1:numParticles
            if ~isempty(Particles{CurrentChannel}(i).Nucleus)
                AssignedNuclei(i)=Particles{CurrentChannel}(i).Nucleus;
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

