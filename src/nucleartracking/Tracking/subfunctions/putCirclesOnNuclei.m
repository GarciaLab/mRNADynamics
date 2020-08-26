function [ ellipse ] = putCirclesOnNuclei(FrameInfo,centers,...
    nFrames, indMitosis, radiusScale, varargin )
%PUTCIRCLESONNUCLEI This function takes up a nuclei structure and the names
% of the images where they are present in order to store a circle with a
% predefined diameter in its ellipse field. This is done to save the time 
% of fitting ellipses by doing it only once all the tracks have been 
% approved.

if nargin == 6
    diameters = varargin{1};
    if numel(diameters) == 1
        diameters = repmat(diameters,nFrames,1);
    end
    if numel(diameters) ~= nFrames
        error('Invalid diameters argument. The diameters vector has to contain as many elements as there are images.')
    end
elseif nargin == 5
    diameters = getDiameters(FrameInfo,nFrames,indMitosis);
else
    error('not enough input arguments');
end

ellipse = cell(nFrames,1);

%(y, x, major axis, minor axis, orientation angle, maxcontourvalue, time,
%particle_id %optionally, schnitz id 
%only the first 4 columns are used, so we'll just set the rest to 0.

for frame = 1:nFrames
    
    nEllipses = size(centers{frame}, 1);
    ellipse{frame} = zeros(nEllipses,8);
    
    for ellipseIndex = 1:size(centers{frame},1)
        
            centroid = fliplr(centers{frame}(ellipseIndex,:));
            majorAxis = 0.5*diameters(frame)*radiusScale;
            minorAxis = 0.5*diameters(frame)*radiusScale;
            orientationAngle = 0;
            maxContourValue = 0;
            time = 0;
            particleID = 0;
            ellipse{frame}(ellipseIndex,:) = [centroid, majorAxis, minorAxis,...
                orientationAngle, maxContourValue, time, particleID];
    end
    
end


end

