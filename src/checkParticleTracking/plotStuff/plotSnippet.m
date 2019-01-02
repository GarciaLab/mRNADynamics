function CurrentSnippet = plotSnippet(snippetFigAxes, rawDataAxes, gaussianAxes, xTrace, ...
    CurrentZIndex, FullSlicePath, Spots, CurrentChannel, CurrentFrame, ...
    CurrentParticleIndex, ExperimentType, intScale, snippet_size, xSize, ... 
    ySize, SnippetEdge, CurrentSnippet)
%PLOTSNIPPET Summary of this function goes here
%   Detailed explanation goes here

if  ~isempty(xTrace) && ~isempty(CurrentZIndex)
    %Get the snippet and the mask, and overlay them
    %(MT, 2018-02-12): lattice data could use this, changed CurrentChannel to coatChannel
    FullSlice=imread(FullSlicePath);
    xSpot = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).xDoG(CurrentZIndex);
    ySpot = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).yDoG(CurrentZIndex);

    if isfield(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex), 'snippet_size') && ~isempty(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).snippet_size)
        snippet_size = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).snippet_size;
        %(MT, 2018-02-12): Hacky fix to get this to run with lattice data -
        %FIX LATER
    elseif strcmpi(ExperimentType,'lattice')
        snippet_size = 13; %pixels
    elseif isfield(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex), 'Snippet')
        try
            snippet_size = floor(size(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).Snippet{1}, 1)/2);
        catch
        end
    end

    CurrentSnippet = double(FullSlice(max(1,ySpot-snippet_size):min(ySize,ySpot+snippet_size),...
        max(1,xSpot-snippet_size):min(xSize,xSpot+snippet_size)));
    imSnippet = mat2gray(CurrentSnippet);
    SnippetEdge=size(CurrentSnippet,1);
    IntegrationRadius = 6*intScale; % this appears to be hard-coded into IdentifySingleSpot
    [xGrid, yGrid] = meshgrid(1:SnippetEdge,1:SnippetEdge);
    rGrid = sqrt((xGrid-ceil(SnippetEdge/2)).^2 + (yGrid-ceil(SnippetEdge/2)).^2);
    SnippetMask = rGrid <= IntegrationRadius;
    IntegrationArea=bwperim(SnippetMask);

    SnippetOverlay=cat(3,IntegrationArea/2 + ...
        +imSnippet,imSnippet,imSnippet);

    imshow(SnippetOverlay,...
        [],'Border','Tight','InitialMagnification',1000, 'Parent', snippetFigAxes)

    hold(snippetFigAxes,'on')

    %this displays the actual snippet for the current z-slice used in
    %intensity calculations. the center may differ from the circle in
    %the overlay figure, which indicates the x-y center of the spot
    %within the brightest z-slice.
    SnippetX=(SnippetEdge-1)/2+1-...
        (Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).xDoG(CurrentZIndex)-...
        Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).xFit(CurrentZIndex));
    SnippetY=(SnippetEdge-1)/2+1-...
        (Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).yDoG(CurrentZIndex)-...
        Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).yFit(CurrentZIndex));
    hold(snippetFigAxes,'off')
else
    imshow(zeros(SnippetEdge), 'Parent', snippetFigAxes)
end

[mesh_y,mesh_x] = meshgrid(1:size(CurrentSnippet,2), 1:size(CurrentSnippet,1));

% Single gaussian function: In future this should be a standalone
% function file to ensure consistency with function used for fitting
singleGaussian = @(params) (params(1).*...
    exp(-(...
    (((cos(params(7)))^2 / (2*params(3)^2) ) + ((sin(params(7)))^2 / 2*params(5)^2))  .* (mesh_x-params(2)).^2 ...
    - 2*((-sin(2*params(7)) / (4*params(3)^2) ) + (sin(2*params(7)) / 4*params(5)^2)) .* (mesh_x-params(2)).*(mesh_y-params(4))...
    + (((sin(params(7)))^2 / (2*params(3)^2) ) + ((cos(params(7)))^2 / 2*params(5)^2)).* (mesh_y-params(4)).^2 ...
    )))...
    + params(6) - CurrentSnippet;

if ~isempty(xTrace) && ~isempty(CurrentZIndex)
    if isfield(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex),'gaussParams')
        gaussParams = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).gaussParams;

        if ~isempty(gaussParams)
            gaussParams= gaussParams{CurrentZIndex};
            try
                gauss = singleGaussian(gaussParams);
            catch
                %not sure in what situation this fails. -AR
                %9/15/2018
                gauss = NaN;
            end
        else
            gauss = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).gaussSpot{CurrentZIndex};
        end

    elseif isfield(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex), 'gaussSpot')
        gauss = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).gaussSpot{CurrentZIndex};
    else
        error('No Gaussian Fit Params or Gauss Snippet Found. Try Re-running segmentSpots')
    end

    if ~isnan(gauss)
        surf(gaussianAxes, gauss + CurrentSnippet);
    end
    title(gaussianAxes,'Gaussian fit')
    zlimit = max(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).CentralIntensity);
    zlim(gaussianAxes,[0, zlimit]);
    surf(rawDataAxes,CurrentSnippet)
    title(rawDataAxes,'Raw data');
    zlim(rawDataAxes,[0, zlimit]);
else
    cla(gaussianAxes, 'reset')
    cla(rawDataAxes, 'reset')
end
end

