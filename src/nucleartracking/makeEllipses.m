function Ellipses = makeEllipses(probabilityMap, probabilityThreshold)


nuclearMask = probabilityMap > probabilityThreshold;
nFrames = size(nuclearMask, 3);

Ellipses = cell(nFrames, 1);

for f = 1:nFrames
    
     nuclearMaskFrame= nuclearMask(:, :, f);
    probabilityMapFrame = probabilityMap(:, :, f);
    
    ellipseStats = regionprops(nuclearMaskFrame, probabilityMapFrame, {...
        'WeightedCentroid',...
        'MajorAxisLength',...
        'MinorAxisLength',...
        'Orientation'});
    
    % Ellipses are definied as
    % x, y, a, b, theta, maxcontourvalue, time, particle_id)
    Ellipses{f, 1} = zeros(length(ellipseStats), 8);   
    
    mjx = [ellipseStats.MajorAxisLength];
    mnx = [ellipseStats.MinorAxisLength];
    ori =  [ellipseStats.Orientation];
    wcs= zeros(length(ellipseStats), 2);
    
    for region = 1:length(ellipseStats)
        wcs(region, :) = ellipseStats(region).WeightedCentroid;
        Ellipses{f, 1}(region, :) = [wcs(region,2),wcs(region,1),mjx(region),mnx(region),ori(region),0,0,0];   
    end
    
    
end

end