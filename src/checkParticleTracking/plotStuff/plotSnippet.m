function [CurrentSnippet, himage] = plotSnippet(snippetFigAxes, rawDataAxes, gaussianAxes, xTrace, ...
    CurrentZIndex, FullSlicePath, Spots, CurrentChannel, CurrentFrame, ...
    CurrentParticleIndex, ExperimentType, intScale, snippet_size, xSize, ... 
    ySize, SnippetEdge, FrameInfo, CurrentSnippet, himage)
%PLOTSNIPPET Summary of this function goes here
%   Detailed explanation goes here

% Spots = castStructNumbersToDoubles(Spots);
intScale = double(intScale);

    if  ~isempty(xTrace) && ~isempty(CurrentZIndex)
        %Get the snippet and the mask, and overlay them
        %(MT, 2018-02-12): lattice data could use this, changed CurrentChannel to coatChannel
        FullSlice=imread(FullSlicePath);
        xSpot = double(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).xDoG(CurrentZIndex));
        ySpot = double(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).yDoG(CurrentZIndex));

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
        IntegrationRadius = 6*intScale; % this appears to be hard-coded into IdentifySingleSpot
        [xGrid, yGrid] = meshgrid(1:SnippetEdge,1:SnippetEdge);
        rGrid = sqrt((xGrid-ceil(SnippetEdge/2)).^2 + (yGrid-ceil(SnippetEdge/2)).^2);
        IntegrationArea= rGrid < IntegrationRadius & (rGrid+1) >= IntegrationRadius;

        SnippetOverlay=cat(3,IntegrationArea/2 + ...
            +imSnippet,imSnippet,imSnippet);

        if ~isempty(himage)
            himage = imshow(SnippetOverlay,...
                [],'Border','Tight','InitialMagnification',1000, 'Parent', snippetFigAxes)
        else
            himage.CData = SnippetOverlay;
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
        SnippetX=(SnippetEdge-1)/2+1-...
            double(((Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).xDoG(CurrentZIndex)))-...
            double(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).xFit(CurrentZIndex)));
        SnippetY=(SnippetEdge-1)/2+1-...
           double( (Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).yDoG(CurrentZIndex))-...
           double( Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).yFit(CurrentZIndex)));
        hold(snippetFigAxes,'off')
    else
        imshow(zeros(SnippetEdge), 'Parent', snippetFigAxes)
    end

    if ~isempty(xTrace) && ~isempty(CurrentZIndex)
        if isfield(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex),'gaussParams')
            gaussParams = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).gaussParams;
     
            if ~isempty(gaussParams)
                gaussParams= double(gaussParams{CurrentZIndex});
                if length(gaussParams) == 7
                    gaussParams = [gaussParams, 0, 0];
                end
                try
                    [y,x] = meshgrid(1:size(CurrentSnippet,2), 1:size(CurrentSnippet,1));
                    g = gaussianForSpot(y, x, CurrentSnippet);
                    gauss = g(gaussParams) + CurrentSnippet;
                catch
                    %not sure in what situation this fails. -AR
                    %9/15/2018
                    gauss = NaN;
                end
            else
                gauss = double(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).gaussSpot{CurrentZIndex});
            end

        elseif isfield(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex), 'gaussSpot')
            gauss = double(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).gaussSpot{CurrentZIndex});
        else
            error('No Gaussian Fit Params or Gauss Snippet Found. Try Re-running segmentSpots')
        end

        if ~isnan(gauss)
            surf(gaussianAxes, gauss);
        end
        title(gaussianAxes,'Gaussian fit')
        zlimit = max(double(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).CentralIntensity));
        zlim(gaussianAxes,[0, zlimit]);
        surf(rawDataAxes,CurrentSnippet)
        title(rawDataAxes,'Raw data');
        zlim(rawDataAxes,[0, zlimit]);
        box(gaussianAxes, 'on')
        box(rawDataAxes, 'on')
        %calibrate the axes
        rawDataAxes.Children.XData = (rawDataAxes.Children.XData - max(rawDataAxes.Children.XData)/2)*xyRes;
        rawDataAxes.Children.YData = (rawDataAxes.Children.YData -max(rawDataAxes.Children.YData)/2)*xyRes;
        gaussianAxes.Children.XData = (gaussianAxes.Children.XData -max(gaussianAxes.Children.XData)/2)*xyRes;
        gaussianAxes.Children.YData = (gaussianAxes.Children.YData -max(gaussianAxes.Children.YData)/2)*xyRes;
        xlabel(rawDataAxes, '\mum'); ylabel(rawDataAxes, '\mum');
        xlabel(gaussianAxes, '\mum'); ylabel(gaussianAxes, '\mum');
    else
        cla(gaussianAxes, 'reset')
        cla(rawDataAxes, 'reset')
    end
end

