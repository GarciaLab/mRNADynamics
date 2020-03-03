function [CurrentSnippet, snipImageHandle] = plotSnippet(snippetFigAxes, rawDataAxes, gaussianAxes, xTrace, ...
    cptState, ExperimentType, snippet_size, xSize, ... 
    ySize, SnippetEdge, CurrentSnippet, snipImageHandle, pixelSize)
%PLOTSNIPPET Summary of this function goes here
%   Detailed explanation goes here

plotSnippet(snippetFigAxes, rawDataAxes, gaussianAxes, xTrace, ...
        cptState.CurrentZIndex, cptState.ImageMat, cptState.Spots, cptState.CurrentChannel, cptState.CurrentFrame, ...
        cptState.CurrentParticleIndex, ExperimentType, snippet_size, xSize, ...
        ySize, SnippetEdge, cptState.FrameInfo, CurrentSnippet, snipImageHandle, pixelSize);

% Spots = castStructNumbersToDoubles(Spots);

scale = 10; %magnification of snippet

    if  ~isempty(xTrace) && ~isempty(cptState.CurrentZIndex)
        % Get the snippet and the mask, and overlay them
        % (MT, 2018-02-12): lattice data could use this, changed CurrentChannel to coatChannel
        xSpot = cptState.getCurrentXDoG();
        ySpot = cptState.getCurrentYDoG();

        if isfield(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex), 'snippet_size') && ~isempty(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).snippet_size)
            
            snippet_size = double(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).snippet_size);
            %(MT, 2018-02-12): Hacky fix to get this to run with lattice data -
            %FIX LATER
        elseif strcmpi(ExperimentType,'lattice')
            snippet_size = 13; %pixels
        end
        snippet_size = snippet_size(1);
        CurrentSnippet = double(FullSlice(max(1,ySpot-snippet_size):min(ySize,ySpot+snippet_size),...
            max(1,xSpot-snippet_size):min(xSize,xSpot+snippet_size)));
        imSnippet = mat2gray(CurrentSnippet);
        SnippetEdge=size(CurrentSnippet,1);
        IntegrationRadius = 6*ceil(sqrt(212/pixelSize)); %integrate 109 pixels around the spot with 212nm pixel size
        [xGrid, yGrid] = meshgrid(1:SnippetEdge,1:SnippetEdge);
        rGrid = sqrt((xGrid-ceil(SnippetEdge/2)).^2 + (yGrid-ceil(SnippetEdge/2)).^2);
        IntegrationArea= rGrid < IntegrationRadius & (rGrid+1) >= IntegrationRadius;

        SnippetOverlay=cat(3,IntegrationArea/2 + ...
            +imSnippet,imSnippet,imSnippet);

        if isempty(snipImageHandle) 
            snipImageHandle = imshow(SnippetOverlay,...
                [],'Border','Tight','InitialMagnification',scale*100, 'Parent', snippetFigAxes);
        else
            snipImageHandle.CData = SnippetOverlay;
        end

        hold(snippetFigAxes,'on')

        scaleBarLength = 6; 
        scaleBarX = 29; scaleBarY = 30;
        plot(snippetFigAxes,[scaleBarX; scaleBarX+scaleBarLength], [scaleBarY;scaleBarY], '-w','LineWidth', 3) %scale bar
        xyRes = FrameInfo(1).PixelSize;
        scaleBarTxt = [num2str(xyRes*(scaleBarLength+1),3), ' \mum'];
        text(snippetFigAxes,scaleBarX-5,scaleBarY+5, scaleBarTxt, 'Color', 'white', 'FontSize', 5)


        %this displays the actual snippet for the current z-slice used in
        %intensity calculations. the center may differ from the circle in
        %the overlay figure, which indicates the x-y center of the spot
        %within the brightest z-slice.

        %JP says: re-check operator precedence here
        SnippetX = (SnippetEdge - 1)/2+1 - cptState.getCurrentXDoG() - cptState.getCurrentXFit();
        SnippetY = (SnippetEdge - 1)/2+1 - cptState.getCurrentYDoG() - cptState.getCurrentYFit();
        hold(snippetFigAxes,'off')
    else
        if isempty(snipImageHandle)
            snipImageHandle = imshow(zeros(SnippetEdge), 'Parent', snippetFigAxes);
        else
            scale = 10;
            snipImageHandle.CData = zeros(SnippetEdge*scale, SnippetEdge*scale, 3);
        end
    end
end
