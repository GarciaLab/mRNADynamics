function [snipFig, snipAxes, snipHandles]  =...
    plotNucleusSnippet(snipFig, snipAxes, snipHandles,...,...
    cntState, hisMat, xSize, ySize)
    hasInputChannels = ~isempty(cntState.liveExperiment.inputChannels);
    hisImage = hisMat(:, :, cntState.CurrentFrame);
    fr_idx = find(cntState.Frames == cntState.CurrentFrame);
    if hasInputChannels
        maxz_idx = cntState.MaxZ(fr_idx);
        medz_idx = cntState.MedZ(fr_idx);
        midmedz_idx = cntState.MidMedZ(fr_idx);
    end
    pixelSize = cntState.FrameInfo(1).PixelSize; % (um)
    nucleusDiameter = getDefaultParameters(cntState.FrameInfo,['d' num2str(cntState.schnitzcells(cntState.CurrentNucleus).cycle)]);
    nucleusDiameter_pixels = nucleusDiameter/pixelSize;
    
    scale = 25; %Corresponds to 50% magnification of the snippet 
    snippet_size = double(round(nucleusDiameter_pixels));
    
  
    xSchnitz = cntState.getCurrentX();
    ySchnitz = cntState.getCurrentY();
    
    if hasInputChannels
        imSnippet = zeros(2*snippet_size + 1, 2*snippet_size + 1, 'double');
        imMedSnippet = zeros(2*snippet_size + 1, 2*snippet_size + 1, 'double');
        imMidMedSnippet = zeros(2*snippet_size + 1, 2*snippet_size + 1, 'double');
    end
    imHisSnippet = zeros(2*snippet_size + 1, 2*snippet_size + 1, 'double');
    if (~isempty(xSchnitz)) & (~isempty(ySchnitz))
        if ySchnitz-snippet_size < 1
            rmin_idx = snippet_size-ySchnitz+2;
        else
            rmin_idx = 1;
        end
        if xSchnitz-snippet_size < 1
            cmin_idx = snippet_size-xSchnitz+2;
        else
            cmin_idx = 1;
        end
        if ySchnitz+snippet_size > ySize
            rmax_idx = rmin_idx + 2*snippet_size - (ySchnitz + snippet_size-ySize);
        else
            rmax_idx = size(imHisSnippet, 1);
        end
        if xSchnitz+snippet_size > xSize
            cmax_idx = cmin_idx + 2*snippet_size - (xSchnitz + snippet_size-xSize);
        else
            cmax_idx =size(imHisSnippet, 2);
        end
        if hasInputChannels
            imMidMedSnippet(rmin_idx:rmax_idx, cmin_idx:cmax_idx) = ...
                mat2gray( double(cntState.ImageMat(max(1,ySchnitz-snippet_size):min(ySize,ySchnitz+snippet_size),...
                        max(1,xSchnitz-snippet_size):min(xSize,xSchnitz+snippet_size))));

            imSnippet(rmin_idx:rmax_idx, cmin_idx:cmax_idx) = ...
                mat2gray( double(cntState.MaxImageMat(max(1,ySchnitz-snippet_size):min(ySize,ySchnitz+snippet_size),...
                        max(1,xSchnitz-snippet_size):min(xSize,xSchnitz+snippet_size))));

            imMedSnippet(rmin_idx:rmax_idx, cmin_idx:cmax_idx) = ...
                mat2gray( double(cntState.MedImageMat(max(1,ySchnitz-snippet_size):min(ySize,ySchnitz+snippet_size),...
                        max(1,xSchnitz-snippet_size):min(xSize,xSchnitz+snippet_size))));
        end

        imHisSnippet(rmin_idx:rmax_idx, cmin_idx:cmax_idx) = ...
            mat2gray( double(hisImage(max(1,ySchnitz-snippet_size):min(ySize,ySchnitz+snippet_size),...
                    max(1,xSchnitz-snippet_size):min(xSize,xSchnitz+snippet_size))));
    end
    IntegrationRadius = 2/pixelSize; % 6*ceil(sqrt(212/pixelSize)); %integrate 109 pixels around the spot with 212nm pixel size
    [xGrid, yGrid] = meshgrid(1:2*snippet_size+1,1:2*snippet_size+1);
    rGrid = sqrt((xGrid-ceil(snippet_size)).^2 + (yGrid-ceil(snippet_size)).^2);
    IntegrationArea= rGrid < IntegrationRadius & (rGrid+1) >= IntegrationRadius;
    
    if hasInputChannels
        SnippetOverlay=cat(3,IntegrationArea/2 + ...
            +imSnippet,imSnippet,imSnippet);
        SnippetMedOverlay=cat(3,IntegrationArea/2 + ...
            +imMedSnippet,imMedSnippet,imMedSnippet);
        SnippetMidMedOverlay=cat(3,IntegrationArea/2 + ...
            +imMidMedSnippet,imMidMedSnippet,imMidMedSnippet);
    end
    HisSnippetOverlay=cat(3,IntegrationArea/2 + ...
        +imHisSnippet,imHisSnippet,imHisSnippet);
    %% 

    snipHandles{1}.CData = HisSnippetOverlay;
    if hasInputChannels
        snipHandles{2}.CData = SnippetMidMedOverlay;
        snipHandles{3}.CData = SnippetOverlay;
        snipHandles{4}.CData = SnippetMedOverlay;
    end
end
