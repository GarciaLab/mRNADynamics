function schnitzcells = setParentNucleus(schnitzcells, ...
    CurrentFrame, CurrentChannel, CurrentParticle, Particles)
%SETPARENTNUCLEUS Summary of this function goes here
%   Detailed explanation goes here

disp('Select the mother nucleus or press enter to delete mother information')
NewNuclei=ginput(1);

if isempty(NewNuclei)
    error('Write the disconnect schnitz part')
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
        Distance=sqrt((NewNuclei(i,1)-xPosSuspect).^2+(NewNuclei(i,2)-yPosSuspect).^2);
        [MinValue,ClosestNucleusIndex]=min(Distance);

        ClickedSchnitz(i)=SchnitzSuspect(ClosestNucleusIndex);
    end

    %Now look at the different cases. Note that I don't have a good
    %way to fix the parent nucleus itself. This might be a bad idea
    %after all



    if length(ClickedSchnitz)==1
        schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).P=ClickedSchnitz;

    elseif isempty(ClickedSchnitz)
        schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).P=0;
    end
end
end

