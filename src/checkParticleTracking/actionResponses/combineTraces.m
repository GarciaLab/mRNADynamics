function [PreviousParticle, Particles] =...
    combineTraces(Spots, CurrentChannelIndex, CurrentFrame, Particles, CurrentParticle)
%COMBINETRACES Summary of this function goes here
%   Detailed explanation goes here

PreviousParticle=0;

exitConnectFlag = 0;

particlesExistInFrame = length(Spots{CurrentChannelIndex}(CurrentFrame).Fits);
currentParticleExistsInCurrentFrame =...
    sum(Particles{CurrentChannelIndex}(CurrentParticle).Frame==CurrentFrame);

if ~currentParticleExistsInCurrentFrame && particlesExistInFrame

    [ConnectPositionx,ConnectPositiony]=ginput(1);
    ConnectPosition = [ConnectPositionx,ConnectPositiony];

    if ~isempty(ConnectPosition)
        % find index of the particle we want to add (a.k.a output particle) to current
        % particle (current particle)
        [ParticleOutput,~]=FindClickedParticle(ConnectPosition,CurrentFrame,...
            Spots{CurrentChannelIndex},Particles{CurrentChannelIndex});

        %Check that the clicked particle doesn't exist in a previous
        %frame, that there is no overlap of frames. If it does
        %exist in a previous frame we will have to disconnect it.

        clickedParticleExistsInAnyPreviousFrame = sum(Particles{CurrentChannelIndex}(ParticleOutput).Frame<CurrentFrame);


        if clickedParticleExistsInAnyPreviousFrame

            msgbox('this button doesn''t currently support adding traces in this direction. try changing to this particle and then adding to the future particle.')
            exitConnectFlag = 1;

        end

        if ~exitConnectFlag
            %Check that there is no overlap. If so, split current particle
            overlap=0;
            for i=1:length(Particles{CurrentChannelIndex}(ParticleOutput).Frame)
                for j=1:length(Particles{CurrentChannelIndex}(CurrentParticle).Frame)
                    if Particles{CurrentChannelIndex}(ParticleOutput).Frame(i)==Particles{CurrentChannelIndex}(CurrentParticle).Frame(j)
                        overlap=1;
                    end
                end
            end

            if overlap
                %Disconnect the clicked particle
                Particles{CurrentChannelIndex}=SeparateParticleTraces(CurrentParticle,CurrentFrame,Particles{CurrentChannelIndex});

                %If the clicked particle has an index larger than that
                %of the current particle we also need to
                %move the index of the clicked particle by one.
                if ParticleOutput>CurrentParticle
                    ParticleOutput=ParticleOutput+1;
                end
            end



            Particles{CurrentChannelIndex}=JoinParticleTraces(CurrentParticle,ParticleOutput,Particles{CurrentChannelIndex});
            %Deals with the indexing changing because of the removal of
            %the old particle.
            if ParticleOutput<CurrentParticle
                CurrentParticle=CurrentParticle-1;
            end
            %Sort the frames within the particle. This is useful if we
            %connected to a particle that came before.
            [SortedFrame,Permutations]=sort(Particles{CurrentChannelIndex}(CurrentParticle).Frame);
            Particles{CurrentChannelIndex}(CurrentParticle).Frame=Particles{CurrentChannelIndex}(CurrentParticle).Frame(Permutations);
            Particles{CurrentChannelIndex}(CurrentParticle).Index=Particles{CurrentChannelIndex}(CurrentParticle).Index(Permutations);
            Particles{CurrentChannelIndex}(CurrentParticle).FrameApproved=Particles{CurrentChannelIndex}(CurrentParticle).FrameApproved(Permutations);

        end
    end

else

    [ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'b', 'linewidth',1);
    ConnectPosition = [ConnectPositionx,ConnectPositiony];

    [ParticleOutput,~]=FindClickedParticle(ConnectPosition,CurrentFrame,Spots{CurrentChannelIndex},Particles{CurrentChannelIndex});

    %If it's an independent particle swap it with the frame in the
    %current particle
    if (length(Particles{CurrentChannelIndex}(ParticleOutput).Frame)==1)&&...
            (sum(Particles{CurrentChannelIndex}(ParticleOutput).Frame==CurrentFrame)==1)

        ParticleTemp=Particles{CurrentChannelIndex}(ParticleOutput);

        %Copy the particle out
        Particles{CurrentChannelIndex}(ParticleOutput).Index=...
            Particles{CurrentChannelIndex}(CurrentParticle).Index(Particles{CurrentChannelIndex}(CurrentParticle).Frame==CurrentFrame);

        %Copy the new particle in
        Particles{CurrentChannelIndex}(CurrentParticle).Index(Particles{CurrentChannelIndex}(CurrentParticle).Frame==CurrentFrame)=...
            ParticleTemp.Index;
        Particles{CurrentChannelIndex}(CurrentParticle).FrameApproved(Particles{CurrentChannelIndex}(CurrentParticle).Frame==CurrentFrame)=1;
    else
        disp('Cannnot connect to two particles!')
    end


end
end

