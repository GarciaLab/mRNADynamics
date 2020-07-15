function ellipseFrame = adjustNuclearContours(ellipseFrame, HisImage, PixelSize_um)
    
    nuclearMask = makeNuclearMask(ellipseFrame,...
        size(HisImage), 1.0);

    %parameters i've found to be broadly applicable
    sigmaSnakes_um = .5;
    mu = .1; %weight of length term for chen vese  algorithm. honestly don't know what this controls
    min_rad_um = .1; % set min and max acceptable area for nucleus segmentation
    max_rad_um = 6; %this needs to be 6um for nc12. 4um for nc14
    nIterSnakes = 100;
    maxAspectRatio = 4;

    minArea_px = round(pi*(min_rad_um ./ PixelSize_um).^2);
    maxArea_px = round(pi*(max_rad_um ./ PixelSize_um).^2);
    sigmaSnakes_px = sigmaSnakes_um / PixelSize_um;
    areaFilter = [minArea_px, maxArea_px];

    nuclearMaskRefined = gather( chenvese( ...
        imgaussfilt( gpuArrayMaybe(HisImage), sigmaSnakes_px),...
        gpuArrayMaybe(nuclearMask), nIterSnakes, mu, 'chan') );

    nuclearMaskSuperRefined = wshed(~nuclearMaskRefined, 'checkPolarity', false);

    [~, ellipseFrame] = fitEllipsesToNuclei(...
        nuclearMaskSuperRefined, 'areaFilter', areaFilter,...
        'maxAspectRatio', maxAspectRatio);

    %validate sizes. the ellipse masker handles
    %very large objects poorly

    if ~isempty(ellipseFrame)
        largeAxisIndex = ellipseFrame(:, 3) > min(size(HisImage))...
            | ellipseFrame(3) > min(size(HisImage));
        ellipseFrame(largeAxisIndex, 3) = median(ellipseFrame(:, 3)); 
        ellipseFrame(largeAxisIndex, 4) = median(ellipseFrame(:, 4));
    end

    ellipseFrame(:, 6:9) = zeros(size(ellipseFrame, 1), 4);    


end