function schnitzcells=CheckMultipleDaughters(schnitzcells,DaughterSchnitz,NewParentSchnitz)

%Checks that the same daughter schnitz is not assigned to multiple parent
%schnitzes.


schnitzcells(NewParentSchnitz)

for i=1:length(schnitzcells)
    if ~(NewParentSchnitz==i)
        if schnitzcells(i).E==DaughterSchnitz
            schnitzcells(i).E=0;
        elseif schnitzcells(i).D==DaughterSchnitz
            schnitzcells(i).D=0;
        end
    end
end

        