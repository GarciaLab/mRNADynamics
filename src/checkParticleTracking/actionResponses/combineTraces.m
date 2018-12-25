function [PreviousParticle, Particles] =...
    combineTraces(Spots, CurrentChannel, CurrentFrame, Particles, CurrentParticle)
%COMBINETRACES Summary of this function goes here
%   Detailed explanation goes here

PreviousParticle=0;

exitConnectFlag = 0;

particlesExistInFrame = length(Spots{CurrentChannel}(CurrentFrame).Fits);
currentParticleExistsInCurrentFrame = sum(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame);

if ~currentParticleExistsInCurrentFrame && particlesExistInFrame

    [ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'b', 'linewidth',1);
    ConnectPosition = [ConnectPositionx,ConnectPositiony];

    if ~isempty(ConnectPosition)
        % find index of the particle we want to add (a.k.a output particle) to current
        % particle (current particle)
        [ParticleOutput,~]=FindClickedParticle(ConnectPosition,CurrentFrame,Spots{CurrentChannel},Particles{CurrentChannel});

        %Check that the clicked particle doesn't exist in a previous
        %frame, that there is no overlap of frames. If it does
        %exist in a previous frame we will have to disconnect it.

        clickedParticleExistsInAnyPreviousFrame = sum(Particles{CurrentChannel}(ParticleOutput).Frame<CurrentFrame);


        if clickedParticleExistsInAnyPreviousFrame

            msgbox('this button doesn''t currently support adding traces in this direction. try changing to this particle and then adding to the future particle.')
            exitConnectFlag = 1;
            %
            %                     %Disconnect the clicked particle
            %                     Particles{CurrentChannel}=SeparateParticleTraces(ParticleOutput,CurrentFrame,Particles{CurrentChannel});
            %                     ParticleOutput=ParticleOutput+1;
            %
            %                     %If the current particle has an index larger than that
            %                     %of the clicked particle (ParticleOutput) we also need to
            %                     %move the index of the current Particle by one.
            %                     if ParticleOutput<CurrentParticle
            %                     	CurrentParticle=CurrentParticle+1;
            %                     end

        end

        if ~exitConnectFlag
            %Check that there is no overlap. If so, split current particle
            overlap=0;
            for i=1:length(Particles{CurrentChannel}(ParticleOutput).Frame)
                for j=1:length(Particles{CurrentChannel}(CurrentParticle).Frame)
                    if Particles{CurrentChannel}(ParticleOutput).Frame(i)==Particles{CurrentChannel}(CurrentParticle).Frame(j)
                        overlap=1;
                    end
                end
            end

            if overlap
                %Disconnect the clicked particle
                Particles{CurrentChannel}=SeparateParticleTraces(CurrentParticle,CurrentFrame,Particles{CurrentChannel});

                %If the clicked particle has an index larger than that
                %of the current particle we also need to
                %move the index of the clicked particle by one.
                if ParticleOutput>CurrentParticle
                    ParticleOutput=ParticleOutput+1;
                end
            end



            Particles{CurrentChannel}=JoinParticleTraces(CurrentParticle,ParticleOutput,Particles{CurrentChannel});
            %Deals with the indexing changing because of the removal of
            %the old particle.
            if ParticleOutput<CurrentParticle
                CurrentParticle=CurrentParticle-1;
            end
            %Sort the frames within the particle. This is useful if we
            %connected to a particle that came before.
            [SortedFrame,Permutations]=sort(Particles{CurrentChannel}(CurrentParticle).Frame);
            Particles{CurrentChannel}(CurrentParticle).Frame=Particles{CurrentChannel}(CurrentParticle).Frame(Permutations);
            Particles{CurrentChannel}(CurrentParticle).Index=Particles{CurrentChannel}(CurrentParticle).Index(Permutations);
            Particles{CurrentChannel}(CurrentParticle).FrameApproved=Particles{CurrentChannel}(CurrentParticle).FrameApproved(Permutations);

        end
    end

else

    [ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'b', 'linewidth',1);
    ConnectPosition = [ConnectPositionx,ConnectPositiony];

    [ParticleOutput,~]=FindClickedParticle(ConnectPosition,CurrentFrame,Spots{CurrentChannel},Particles{CurrentChannel});

    %If it's an independent particle swap it with the frame in the
    %current particle
    if (length(Particles{CurrentChannel}(ParticleOutput).Frame)==1)&&...
            (sum(Particles{CurrentChannel}(ParticleOutput).Frame==CurrentFrame)==1)

        ParticleTemp=Particles{CurrentChannel}(ParticleOutput);

        %Copy the particle out
        Particles{CurrentChannel}(ParticleOutput).Index=...
            Particles{CurrentChannel}(CurrentParticle).Index(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame);

        %Copy the new particle in
        Particles{CurrentChannel}(CurrentParticle).Index(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)=...
            ParticleTemp.Index;
        Particles{CurrentChannel}(CurrentParticle).FrameApproved(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)=1;
    else
        disp('Cannnot connect to two particles!')
    end


end
end

