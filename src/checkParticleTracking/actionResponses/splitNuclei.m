function [Particles, PreviousParticle, schnitzcells] = splitNuclei(schnitzcells, ...
    CurrentFrame, CurrentChannel, CurrentParticle, Particles)
%SPLITNUCLEI Summary of this function goes here
%   Detailed explanation goes here

PreviousParticle=0;
        
disp('Select one/two daughter nuclei and press ENTER or just press ENTER to terminate lineage')
NewNuclei=ginput(2);

if isempty(NewNuclei)
    warning('Write the disconnect schnitz part')
else
    [NClicks,~]=size(NewNuclei);

    ClickedSchnitz=[];
    %Find the nuclei/schnitz that we clicked on
    for i=1:NClicks
        %Find which schnitz this corresponds to
        SchnitzSuspect=[];
        xPosSuspect=[];
        yPosSuspect=[];
        for j=1:length(schnitzcells)
            if sum(schnitzcells(j).frames==CurrentFrame)
                SchnitzSuspect=[SchnitzSuspect,j];
                if (~isempty(schnitzcells(j).cenx(find((schnitzcells(j).frames)==CurrentFrame))))&...
                        (~isempty(schnitzcells(j).ceny(find((schnitzcells(j).frames)==CurrentFrame))))
                    xPosSuspect=[xPosSuspect,...
                        schnitzcells(j).cenx(find((schnitzcells(j).frames)==CurrentFrame))];
                    yPosSuspect=[yPosSuspect,...
                        schnitzcells(j).ceny(find((schnitzcells(j).frames)==CurrentFrame))];
                else
                    xPosSuspect=[xPosSuspect,inf];
                    yPosSuspect=[yPosSuspect,inf];
                end
            end
        end

        %Find the closest one to the point where we clicked
        Distance=sqrt(double((NewNuclei(i,1)-xPosSuspect).^2+(NewNuclei(i,2)-yPosSuspect).^2));
        [MinValue,ClosestNucleusIndex]=min(Distance);

        ClickedSchnitz(i)=SchnitzSuspect(ClosestNucleusIndex);
    end

    %Now look at the different cases


    %       Click on one nucleus + ENTER: Continue the schnitz with that nucleus.
    %       Click on two nuclei: Split the current nucleus into two daughter
    %       nuclei.
    %       Click on the same nucleus twice: Split the current nucleus, but
    %       with only one daughter nucleus.


    if length(ClickedSchnitz)==1
        if Particles{CurrentChannel}(CurrentParticle).Nucleus==ClickedSchnitz %Split the lineage
            [Particles{CurrentChannel},schnitzcells]=SplitSchnitz(Particles{CurrentChannel},schnitzcells,...
                CurrentFrame,...
                CurrentParticle);
        else
            try
                [Particles{CurrentChannel},schnitzcells]=...
                    JoinSchnitz(Particles{CurrentChannel},schnitzcells,Particles{CurrentChannel}(CurrentParticle).Nucleus,...
                    ClickedSchnitz,CurrentFrame);
            catch
                disp('Error in JoinSchnitz')
            end
        end
    elseif length(ClickedSchnitz)==2
        if ClickedSchnitz(1)~=ClickedSchnitz(2)
            [Particles{CurrentChannel},schnitzcells]=SplitSchnitzDaughters(Particles{CurrentChannel},schnitzcells,...
                CurrentFrame,...
                Particles{CurrentChannel}(CurrentParticle).Nucleus,ClickedSchnitz(1),ClickedSchnitz(2));
        else
            [Particles{CurrentChannel},schnitzcells]=SplitSchnitzDaughters(Particles{CurrentChannel},schnitzcells,...
                CurrentFrame,...
                Particles{CurrentChannel}(CurrentParticle).Nucleus,ClickedSchnitz(1),0);
        end
    else
        error('Too many cells selected')
    end
end

%Check if there are any issues with cellno getting lost in the
%schnitz
for SchnitzN=1:length(schnitzcells)
    if length(schnitzcells(SchnitzN).frames)~=length(schnitzcells(SchnitzN).cellno)
        warning(['Problem with schnitz ',num2str(SchnitzN)])
        keyboard
    end
end
end

