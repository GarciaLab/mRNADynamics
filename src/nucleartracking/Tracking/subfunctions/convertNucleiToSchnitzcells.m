function [ schnitzcells, varargout ] =...
    convertNucleiToSchnitzcells( nuclei, approvedCenters,...
    approvedSchnitzcell, previousSchnitzcell)
%CONVERT MAPPING AND CENTER TO SCHNITZCELLS 

approvedSchnitz = false(numel(nuclei),1);


for nucleus = 1:numel(nuclei)

    schnitzcells(nucleus).P = uint16(nuclei(nucleus).P);
    schnitzcells(nucleus).E = uint16(nuclei(nucleus).E);
    schnitzcells(nucleus).D = uint16(nuclei(nucleus).D);
    schnitzcells(nucleus).frames = uint16(find(nuclei(nucleus).indXY > 0));
    for jj = 1:numel(schnitzcells(nucleus).frames)
        frameInd = schnitzcells(nucleus).frames(jj);
        schnitzcells(nucleus).cenx(jj) = uint16(nuclei(nucleus).position(frameInd,2));
        schnitzcells(nucleus).ceny(jj) = uint16(nuclei(nucleus).position(frameInd,1));
        schnitzcells(nucleus).len(jj) = NaN;
        schnitzcells(nucleus).cellno(jj) = uint16(nuclei(nucleus).indXY(frameInd));
        
        if exist('approvedCenters','var') &&...
                approvedCenters{frameInd}(nuclei(nucleus).indXY(frameInd))
            approvedSchnitz(nucleus) = true;
        else
            approvedSchnitz(nucleus) = false;
        end
    end
        
end
varargout{1} = approvedSchnitz;

if exist('approvedSchnitzcell','var')
    map = [];
    for nucleus = 1:numel(previousSchnitzcell)
        if approvedSchnitzcell(nucleus)
            for jj = 1:numel(schnitzcells)
                if schnitzcells(jj).frames(1) == previousSchnitzcell(nucleus).frames(1)...
                        && schnitzcells(jj).cenx(1) == previousSchnitzcell(nucleus).cenx(1)...
                        && schnitzcells(jj).ceny(1) == previousSchnitzcell(nucleus).ceny(1)
                    
                    map = [map; nucleus jj];
                    
                    break
                end
            end
        end
    end
    
    varargout{2} = map;
    
end


end