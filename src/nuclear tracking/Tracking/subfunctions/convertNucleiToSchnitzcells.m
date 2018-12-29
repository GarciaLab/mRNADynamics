function [ schnitzcells, varargout ] = convertNucleiToSchnitzcells( nuclei, approvedCenters, approvedSchnitzcell, previousSchnitzcell)
%CONVERTMAPPINGANDCENTERTOSCHNITZCELLS Summary of this function goes here
%   Detailed explanation goes here

approvedSchnitz = false(numel(nuclei),1);


for j = 1:numel(nuclei)

    schnitzcells(j).P = nuclei(j).P;
    schnitzcells(j).E = nuclei(j).E;
    schnitzcells(j).D = nuclei(j).D;
    schnitzcells(j).frames = find(nuclei(j).indXY > 0);
    for jj = 1:numel(schnitzcells(j).frames)
        frameInd = schnitzcells(j).frames(jj);
        schnitzcells(j).cenx(jj) = nuclei(j).position(frameInd,2);
        schnitzcells(j).ceny(jj) = nuclei(j).position(frameInd,1);
        schnitzcells(j).len(jj) = NaN;
        schnitzcells(j).cellno(jj) = nuclei(j).indXY(frameInd);
        
        if exist('approvedCenters','var') && approvedCenters{frameInd}(nuclei(j).indXY(frameInd))
            approvedSchnitz(j) = true;
        else
            approvedSchnitz(j) = false;
        end
    end
        
end
varargout{1} = approvedSchnitz;

if exist('approvedSchnitzcell','var')
    map = [];
    for j = 1:numel(previousSchnitzcell)
        if approvedSchnitzcell(j)
            for jj = 1:numel(schnitzcells)
                if schnitzcells(jj).frames(1) == previousSchnitzcell(j).frames(1) && schnitzcells(jj).cenx(1) == previousSchnitzcell(j).cenx(1) && schnitzcells(jj).ceny(1) == previousSchnitzcell(j).ceny(1)
                    map = [map; j jj];
                    break
                end
            end
        end
    end
    varargout{2} = map;
end

