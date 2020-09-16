function [PreviousNucleus, schnitzcells] =...
    combineNuclearTraces(CurrentFrame, schnitzcells, Ellipses, CurrentNucleus, FrameInfo, ncFrames)
%COMBINETRACES Summary of this function goes here
%   Detailed explanation goes here

PreviousNucleus=0;

exitConnectFlag = 0;

nucleiExistInFrame = size(Ellipses{CurrentFrame}, 1);
currentNucleusExistsInCurrentFrame =...
    sum(schnitzcells(CurrentNucleus).frames==CurrentFrame);

if ~currentNucleusExistsInCurrentFrame && nucleiExistInFrame

    [ConnectPositionx,ConnectPositiony]=ginput(1);
    ConnectPosition = [ConnectPositionx,ConnectPositiony];

    if ~isempty(ConnectPosition)
        % find index of the particle we want to add (a.k.a output particle) to current
        % particle (current particle)
        [NucleusOutput,~]=FindClickedNucleus(ConnectPosition,CurrentFrame,...
            schnitzcells);
        %Check that the clicked particle doesn't exist in a previous
        %frame, that there is no overlap of frames. If it does
        %exist in a previous frame we will have to disconnect it.

        clickedNucleusExistsInAnyPreviousFrame = sum(schnitzcells(NucleusOutput).frames<CurrentFrame);


        if clickedNucleusExistsInAnyPreviousFrame

            msgbox('this button doesn''t currently support adding traces in this direction. try changing to this particle and then adding to the future particle.')
            exitConnectFlag = 1;

        end

        if ~exitConnectFlag
            %Check that there is no overlap. If so, split current particle
            overlap=0;
            for i=1:length(schnitzcells(NucleusOutput).frames)
                for j=1:length(schnitzcells(CurrentNucleus).frames)
                    if schnitzcells(NucleusOutput).frames(i)==schnitzcells(CurrentNucleus).frames(j)
                        overlap=1;
                    end
                end
            end

            if overlap
                %Disconnect the clicked particle
                schnitzcells=SeparateNucleusTraces(CurrentNucleus,CurrentFrame,schnitzcells,...
                    FrameInfo, ncFrames);

                %If the clicked particle has an index larger than that
                %of the current particle we also need to
                %move the index of the clicked particle by one.
                if NucleusOutput>CurrentNucleus
                    NucleusOutput=NucleusOutput+1;
                end
            end



            schnitzcells=JoinNuclearTraces(CurrentNucleus,NucleusOutput,schnitzcells, FrameInfo, ncFrames);
            %Deals with the indexing changing because of the removal of
            %the old particle.
            if NucleusOutput<CurrentNucleus
                CurrentNucleus=CurrentNucleus-1;
            end

        end
    end

else
% 
%     [ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'b', 'linewidth',1);
%     ConnectPosition = [ConnectPositionx,ConnectPositiony];
% 
%     [NucleusOutput,~]=FindClickedNucleus(ConnectPosition,CurrentFrame,...
%             schnitzcells);
%     %If it's an independent particle swap it with the frame in the
%     %current particle
%     if (length(schnitzcells(NucleusOutput).frames)==1)&&...
%             (sum(schnitzcells(NucleusOutput).frames==CurrentFrame)==1)
% 
%         NucleusTemp=schnitzcells(NucleusOutput);
% 
%         %Copy the particle out
%         schnitzcells(NucleusOutput).cellno=...
%             schnitzcells(CurrentNucleus).cellno(schnitzcells(CurrentNucleus).frames==CurrentFrame);
% 
%         %Copy the new particle in
%         schnitzcells(CurrentNucleus).cellno(schnitzcells(CurrentNucleus).frames==CurrentFrame)=...
%             NucleusTemp.cellno;
%         schnitzcells(CurrentNucleus).FrameApproved(schnitzcells(CurrentNucleus).frames==CurrentFrame)=1;
% %     else
    disp('Cannnot connect to two schnitz cells!')
%     end


end
end

