function [ centers, mapping, approved_centers, varargout ] = convertSchnitzcellsToCentersAndMapping( schnitzcells, nFrames, varargin )
%CONVERTSCHNITZTOCENTERSANDMAPPING The function takes a structure of
% schnitzcells and extracts the position and mapping information to rerun
% the tracking.
%
%   nFrames is the total number of frames in the movie.
%                                           
%   The optionnal parameter is a vector of the size of the schnitzcells 
%   structure that corresponds to the 'Approved' state of each schnitzcell.
%   Any value above 0 is considered as approved.

centers = cell(nFrames,1);
mapping = cell(nFrames-1,1);

if nargin > 2
    approvedSchnitzcells = varargin{1};
    approved_centers = cell(size(centers));
else
    approvedSchnitzcells = false(size(schnitzcells));
    approved_centers = cell(size(centers));
end

if nargout > 3
    convert_to_nuclei = true;
    nuclei = struct('P',zeros(660,0),'E',zeros(660,0),'D',zeros(660,0));
else
    convert_to_nuclei = false;
end

for j = 1:numel(schnitzcells)
    schnitz = schnitzcells(j);
    approvedFlag = approvedSchnitzcells(j);
    for jj = 1:numel(schnitz.frames)
        frameInd = schnitz.frames(jj);
        ind_centers = schnitz.cellno(jj);
        % Copy the position of the nucleus
        centers{frameInd}(ind_centers,:) = [schnitz.ceny(jj) schnitz.cenx(jj)];
        % Copy the tracking information only in cases where the schnitzcell was
        % approved.
        if approvedFlag > 0
            approved_centers{frameInd}(ind_centers) = true;
            if jj < numel(schnitz.frames) % Before division, just track the nucleus
                mapping{frameInd}(ind_centers,:) = [schnitz.cellno(jj+1) 0];
            else % During division, track the two daughters.
                if ~isempty(schnitz.E) && schnitz.E ~= 0
                    daughter1 = schnitzcells(schnitz.E).cellno(1);
                else
                    daughter1 = 0;
                end
                
                if ~isempty(schnitz.D) && schnitz.D ~= 0
                    daughter2 = schnitzcells(schnitz.D).cellno(1);
                else
                    daughter2 = 0;
                end
                daughters = [daughter1 daughter2];
                mapping{frameInd}(ind_centers,:) = sort(daughters,'descend');   
            end
        else
            mapping{frameInd}(ind_centers,:) = [0 0];
            approved_centers{frameInd}(ind_centers) = false;
        end
    end
    
end

if convert_to_nuclei
    for j = 1:numel(schnitzcells)
        
        nuclei(j).P = schnitzcells(j).P;
        nuclei(j).E = schnitzcells(j).E;
        nuclei(j).D = schnitzcells(j).D;
        nuclei(j).indXY = nan(nFrames,1);
        nuclei(j).indXY(schnitzcells(j).frames) = schnitzcells(j).cellno;
        nuclei(j).position = nan(nFrames,2);
        nuclei(j).position(schnitzcells(j).frames,:) = [schnitzcells(j).ceny' schnitzcells(j).cenx'];
        
    end
    varargout{1} = nuclei;
end
end

