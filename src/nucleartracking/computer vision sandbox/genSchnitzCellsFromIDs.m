    function schnitzcells = genSchnitzCellsFromIDs(...
        schnitzcells, ids, frameIndex, cenxs, cenys, len)

 
% This is a subfunction for tracknuclei_computervision 

%AR 4/2020
        for k = 1:length(ids)
            if length(schnitzcells) < k
                
                schnitzcells(k).cenx(1) = cenxs(k);
                schnitzcells(k).ceny(1) = cenys(k);
                schnitzcells(k).frames(1) = frameIndex;
                schnitzcells(k).len(1) = len(k);
                
                
            else
                
                schnitzcells(k).cenx(end+1) = cenxs(k);
                schnitzcells(k).ceny(end+1) = cenys(k);
                schnitzcells(k).frames(end+1) = frameIndex;
                schnitzcells(k).len(end+1) = len(k);
                
                
            end
        end
    end