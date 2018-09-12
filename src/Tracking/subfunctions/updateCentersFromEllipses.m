function [ centers ] = updateCentersFromEllipses(Prefix, Ellipses, varargin )
%UPDATECENTERSFROMELLIPSES Takes an Ellipses structure and outputs an 
% updated of the centers structure.
%
% There are two ways to use it :
%
%   centers = updateCentersFromEllipses(Ellipses)
%
%   centers is a new structure that contains only the center coordinates of
%   all ellipses.
%
%   centers = updateCentersFromEllipses(Ellipses,centers)
%
%   updates the centers structure by just adding or removing the centers
%   that were not present in both structures, Ellipses and centers.


space_resolution = getDefaultParameters(Prefix,'space resolution');

THRESHOLD_DISTANCE = 2; % in micrometers. Only used when the old centers strucutre is provided too. Maximum offset allowed between a nucleus's position in the Ellipses struct and the centers struct for them to be considered the same nucleus.
THRESHOLD_DISTANCE = THRESHOLD_DISTANCE/space_resolution;

if nargin > 3
    old_centers_provided = true;
    old_centers = varargin{1};
else
    old_centers_provided = false;
end

nFrames = numel(Ellipses);
centers = cell(size(Ellipses));

if old_centers_provided
    
    % If old_centers are provided, the strategy is to compare old_centers
    % and Ellipses and see which nuclei appeared/disappeared and then copy
    % the ones that are in common and add the ones that appeared.
    for j = 1:nFrames
        if size(old_centers{j},2) < 2
            old_centers{j} = [inf inf];
        end
            dist = pdist2(Ellipses{j}(:,1:2),old_centers{j});
            % Build correspondences
            
            % First attribute nuclei that are mutually the closest
            [mc,ic] = min(dist,[],1);
            [mr,ir] = min(dist,[],2);
            
            indc = sub2ind(size(dist),ic,1:numel(ic));
            indr = sub2ind(size(dist),1:numel(ir),ir');
            
            indc(mc>THRESHOLD_DISTANCE) = [];
            %indr(mr>THRESHOLD_DISTANCE) = [];
            
            Mc = zeros(size(dist));
            Mc(indc) = 1;
            Mr = zeros(size(dist));
            Mr(indr(mr'<THRESHOLD_DISTANCE)) = 1;
            
            ind = find(Mc & Mr);
            [rowInd,colInd] = ind2sub(size(dist),ind);
            
            centers{j} = old_centers{j}(colInd,:);
            
            % Add the ones that were added
            ind_to_add = 1:size(Ellipses{j},1);
            ind_to_add(rowInd) = [];
            for jj = 1:numel(ind_to_add)
                centers{j} = [centers{j}(1:ind_to_add(jj)-1,:); Ellipses{j}(ind_to_add(jj),1:2); centers{j}(ind_to_add(jj):end,:)];
            end
            
    end
else
    
    % Just copy the information from Ellipses
    % I should make this and the rest of the code be able to handle empty
    % frames!
    for j = 1:nFrames
        %if ~isempty(Ellipses{j})
            centers{j} = Ellipses{j}(:,[2,1]);
        %else
        %    centers{j}=[];
        %end
    end
    
    
    
end

