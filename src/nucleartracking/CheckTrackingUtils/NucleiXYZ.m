function [x,y,z, cellnos]=NucleiXYZ(schnitzcells, CurrentFrame)


%Return the X and Y coordinate of the brightest Z of each spot in the
%Spots structure
x = [];
y = [];
z = [];
cellnos = [];
idx = 1;
if ~isempty(schnitzcells)
    for sc=1:length(schnitzcells)
        if ismember(CurrentFrame, schnitzcells(sc).frames)
            x(idx)=double(schnitzcells(sc).cenx(schnitzcells(sc).frames == CurrentFrame));
            y(idx)=double(schnitzcells(sc).ceny(schnitzcells(sc).frames == CurrentFrame));
            if isfield(schnitzcells, 'Fluo')
                if ~isempty(find(schnitzcells(sc).Fluo(schnitzcells(sc).frames == CurrentFrame,:) == ...
                    max(schnitzcells(sc).Fluo(schnitzcells(sc).frames == CurrentFrame,:)), 1))
                    z(idx)=find(schnitzcells(sc).Fluo(schnitzcells(sc).frames == CurrentFrame,:) == ...
                        max(schnitzcells(sc).Fluo(schnitzcells(sc).frames == CurrentFrame,:)), 1);
                else
                    z(idx) = 1;
                end
            end
            cellnos(idx) = schnitzcells(sc).cellno(schnitzcells(sc).frames == CurrentFrame);
            idx = idx + 1;
        end

    end
else
    x=[];
    y=[];
    z=[];
    cellnos = [];
end

