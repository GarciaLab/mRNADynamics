function [ ellipse ] = putCirclesOnNuclei(FrameInfo,centers, nFrames, indMitosis, varargin )
%PUTCIRCLESONNUCLEI This function takes up a nuclei structure and the names
% of the images where they are present in order to store a circle with a
% predefined diameter in its ellipse field. This is done to save the time 
% of fitting ellipses by doing it only once all the tracks have been 
% approved.

if nargin > 4
    diameters = varargin{1};
    if numel(diameters) == 1
        diameters = repmat(diameters,nFrames,1);
    end
    if numel(diameters) ~= nFrames
        error('Invalid diameters argument. The diameters vector has to contain as many elements as there are images.')
    end
else
    diameters = getDiameters(FrameInfo,nFrames,indMitosis);
end

ellipse = cell(nFrames,1);
        
for j = 1:nFrames
    ellipse{j} = zeros(size(centers{j},1),8);
    for jj = 1:size(centers{j},1)
            ellipse{j}(jj,:) = [fliplr(centers{j}(jj,:)) 0.5*diameters(j)*ones(1,2) 0 0 0 0];
    end
end


end

