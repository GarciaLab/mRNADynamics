function Ellipses = makeEllipses(pMap, thresh)


nuclearMask = pMap > thresh;
nFrames = size(nuclearMask, 1);

Ellipses = cell(nFrames, 1);

for f = 1:nFrames
    
    mask = squeeze(nuclearMask(f, :, :)) > thresh;
    pim = squeeze(pMap(f, :, :));
    
    ellipseStats = regionprops(mask, pim, {...
        'WeightedCentroid',...
        'MajorAxisLength',...
        'MinorAxisLength',...
        'Orientation'});
    
    Ellipses{f, 1} = zeros(length(ellipseStats), 8);   %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
    
    mjx = [ellipseStats.MajorAxisLength];
    mnx = [ellipseStats.MinorAxisLength];
    ori =  [ellipseStats.Orientation];
    wcs= zeros(length(ellipseStats), 2);
    
    for r = 1:length(ellipseStats)
        wcs(r, :) = ellipseStats(r).WeightedCentroid;
        Ellipses{f, 1}(r, :) = [wcs(r,2),wcs(r,1),mjx(r),mnx(r),ori(r),0,0,0];     %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
    end
    
    
end

end