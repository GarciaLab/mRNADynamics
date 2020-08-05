function ellipseFrame = adjustNuclearContours(ellipseFrame, HisImage, PixelSize_um)
    
%%
    ellipseFrameEdges = [];
    ellipseFrameInner = []; 
    
    mask = HisImage;
    xDim = size(mask, 2);
    yDim = size(mask, 1);

    
    %if an object is within borderThresh px of the image edge,
    %let's not fit a circle to it and leave it be
    borderThreshold = median( .5*(ellipseFrame(:, 3) + ellipseFrame(:, 4)) );
    border = borderImage(mask);
    borderDist = bwdist(border);

    for k = 1:size(ellipseFrame, 1)
        xSub = max(min(round(abs(ellipseFrame(k, 1))), xDim), 1);
        ySub = max(min(round(abs(ellipseFrame(k, 2))), yDim), 1);

        isNearBorder =  borderDist(ySub, xSub) < borderThreshold;

        if isNearBorder
            ellipseFrameEdges = [ellipseFrameEdges; ellipseFrame(k, :)];
        else
            ellipseFrameInner = [ellipseFrameInner; ellipseFrame(k, :)];
        end
    end
    
    %Add a ninth column so we can append arrays properly. this is where the
    %schnitz correspondence will be located
    ellipseFrameEdges(:, 9) = zeros(size(ellipseFrameEdges, 1), 1);

    
%%    
    nuclearMask = makeNuclearMask(ellipseFrameInner,...
        size(HisImage), 1.0);

    %parameters i've found to be broadly applicable
    sigmaSnakes_um = .5;
    mu = .1; %weight of length term for chen vese  algorithm. honestly don't know what this controls
    min_rad_um = 2; % set min and max acceptable area for nucleus segmentation
    max_rad_um = 8; %this needs to be 6um for nc12. 4um for nc14
    nIterSnakes = 100;
    maxAspectRatio = 4;

    minArea_px = round(pi*(min_rad_um ./ PixelSize_um).^2);
    maxArea_px = round(pi*(max_rad_um ./ PixelSize_um).^2);
    sigmaSnakes_px = sigmaSnakes_um / PixelSize_um;
    areaFilter = [minArea_px, maxArea_px];

%     s1 = 6;
%     filteredHis = gpuArrayWrapper(-filterImage(HisImage,...
%         'Maximum', {s1}));
    filteredHis = imgaussfilt( gpuArrayWrapper(HisImage),...
    sigmaSnakes_px);

    nuclearMaskRefined = gather( chenvese( ...
        filteredHis,...
        gpuArrayWrapper(nuclearMask), nIterSnakes, mu, 'chan') );

    nuclearMaskSuperRefined = wshed(~nuclearMaskRefined, 'checkPolarity', false);

    [~, ellipseFrame] = fitEllipsesToNuclei(...
        nuclearMaskSuperRefined, 'areaFilter', areaFilter,...
        'maxAspectRatio', maxAspectRatio);
    
    ellipseFrame(:, 6:9) = zeros(size(ellipseFrame, 1), 4);  

    ellipseFrame = [ellipseFrame; ellipseFrameEdges];

    %filter nuclei out by size

    if ~isempty(ellipseFrame)
        
        largeAxisIndex = ellipseFrame(:, 3) > max_rad_um/PixelSize_um...
            | ellipseFrame(:, 4) > max_rad_um/PixelSize_um;
        ellipseFrame(largeAxisIndex, :) = [];
        
        smallAxisIndex = ellipseFrame(:, 3) < min_rad_um/PixelSize_um |...
            ellipseFrame(:, 4) < min_rad_um/PixelSize_um;
        ellipseFrame(smallAxisIndex, :) = [];
        
    end
    

end