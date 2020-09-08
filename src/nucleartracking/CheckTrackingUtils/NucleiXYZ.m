function [x,y,z]=NucleiXYZ(schnitzcells, CurrentFrame)


%Return the X and Y coordinate of the brightest Z of each spot in the
%Spots structure

if ~isempty(schnitzcells)
    for sc=1:length(schnitzcells)
        if ismember(CurrentFrames, schnitzcells(sc).frames)
            x(sc)=double(schnitzcells.cenx(schnitzcells(sc).frames == CurrentFrame));
            y(sc)=double(schnitzcells.ceny(schnitzcells(sc).frames == CurrentFrame));
            z(sc)=find(schnitzcells.Fluo(schnitzcells(sc).frames == CurrentFrame,:) == ...
                max(schnitzcells.Fluo(schnitzcells(sc).frames == CurrentFrame,:)));
        end

    end
else
    x=[];
    y=[];
    z=[];
end

