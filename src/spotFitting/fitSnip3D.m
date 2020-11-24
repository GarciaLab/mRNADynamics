function spotsFrame= fitSnip3D(spotsFrame, spotChannel,...
    spotIndex, frame, Prefix,~, ~, nSpots, imStack, displayFigures)

%%

if ischar(Prefix)
    liveExperiment = LiveExperiment(Prefix);
else
    liveExperiment = Prefix;
end

FrameInfo = getFrameInfo(liveExperiment);
preFolder = liveExperiment.preFolder;

if nargin < 9
    movieMat = getMovieMat(liveExperiment);
    if ~isempty(movieMat)
        imStack = movieMat(:, :, :, frame, spotChannel);
    else
        imStack = getMovieFrame(liveExperiment, frame, spotChannel);
    end
end

% extract basic fit parameters
spot = spotsFrame.Fits(spotIndex);
xSize = FrameInfo(1).PixelsPerLine;
ySize = FrameInfo(1).LinesPerFrame;
pixelSize_nm = FrameInfo(1).PixelSize*1000; %nm
zStep_nm = FrameInfo(1).ZStep*1000; %nm
zMax = FrameInfo(1).NumberSlices+2;
snipDepth = uint8(ceil(2500/zStep_nm));

% NL: Need to make this independent of 2D fit info
brightestZPlane = spot.brightestZ;
xSpot = spot.xDoG(spot.z==brightestZPlane);
ySpot = spot.yDoG(spot.z==brightestZPlane);

if isfield(spot, 'snippet_size') && ~isempty(spot.snippet_size)
    snippet_size = spot.snippet_size;
else
    snippet_size = round(1500/pixelSize_nm); % (in pixels)set to be around 1.5 um
end
snippet_size = uint16(snippet_size(1));


zBot = max([1,brightestZPlane - snipDepth]);
zTop = min([zMax, brightestZPlane + snipDepth]);
zRange = zBot:zTop;
xRange = max([1,xSpot-snippet_size]):min([xSize,xSpot+snippet_size]);
yRange = max([1,ySpot-snippet_size]):min([ySize,ySpot+snippet_size]);
snip3D = NaN(numel(yRange),numel(xRange),numel(zBot:zTop));

%%

if exist('imStack', 'var') && ~isempty(imStack)
    
    snip3D = imStack(yRange, xRange, :);
    snip3D = snip3D(:, :, zRange);
    
else
    
    n = 1;
    for z = zRange
        
        FullSlice=imread([preFolder, filesep,liveExperiment.Prefix,'_',iIndex(frame,3)...
            ,'_z' iIndex(z,2) '_ch' iIndex(spotChannel,2) '.tif']);
        
        snip3D(:,:,n) = double(FullSlice(yRange,xRange));
        n = n + 1;
        
    end
    
end

if size(snip3D, 3) == 1
    warning('Snippet is only 2D. Fit will perform very poorly.')
end

%%
xMin = single(min(xRange));
yMin = single(min(yRange));
zMin = single(min(zRange));
% disable two-spot fits for now
if nSpots == 2
    warning('two spot fitting not currently supported. Performing single spot fit')
    nSpots = 1;
end
if nSpots == 2
    [GaussParams1, GaussParams2, offset, GaussIntegralVector, centroid_mean,...
        GaussSE1, GaussSE2, offsetSE, GaussIntegralSEVector, centroid_se] = ...
        fit3DGaussian2spot(snip3D,pixelSize_nm);
    
    % spot 1 position
    spotsFrame.Fits(spotIndex).Spot1Fits3D = single(GaussParams1);
    spotsFrame.Fits(spotIndex).Spot1FitsSE3D = single(GaussSE1);
    x1 = single(GaussParams1(3) + xMin - 1);
    y1 = single(GaussParams1(2) + yMin - 1);
    z1 = single(GaussParams1(4) + zMin - 1);
    spotsFrame.Fits(spotIndex).Spot1Pos3D = single([x1,y1,z1]);
    spotsFrame.Fits(spotIndex).Spot1Int3D = GaussIntegralVector(1);
    spotsFrame.Fits(spotIndex).Spot1IntSE3D = GaussIntegralSEVector(1);
    
    % spot 2 position
    spotsFrame.Fits(spotIndex).Spot2Fits3D = single(GaussParams2);
    spotsFrame.Fits(spotIndex).Spot2CI3D = single(GaussSE2);
    x2 = single(GaussParams2(3) + xMin - 1);
    y2 = single(GaussParams2(2) + yMin - 1);
    z2 = single(GaussParams2(4) + zMin - 1);
    spotsFrame.Fits(spotIndex).Spot2Pos3D = single([x2,y2,z2]);
    spotsFrame.Fits(spotIndex).Spot2Int3D = GaussIntegralVector(2);
    spotsFrame.Fits(spotIndex).Spot2IntSE3D = GaussIntegralSEVector(2);
    
    % combined metrics
    spotsFrame.Fits(spotIndex).gauss3DIntensity = GaussIntegralVector(3);
    spotsFrame.Fits(spotIndex).gauss3DIntensitySE = GaussIntegralSEVector(3);
    spotsFrame.Fits(spotIndex).Offset3D = offset;
    spotsFrame.Fits(spotIndex).Offset3DSE = offsetSE;
    spotsFrame.Fits(spotIndex).GaussPos3D = centroid_mean;
    spotsFrame.Fits(spotIndex).GaussPos3DSE = centroid_se;
    
    % single spot fit
elseif nSpots == 1

    [GaussFit, FitDeltas, GaussIntegral, GaussIntegralSE,...
        GaussIntegralRaw,isSymmetricPositiveDefinite,integrationDims]  = ...
        fit3DGaussianRho(snip3D,[pixelSize_nm zStep_nm]);
    
    spotsFrame.Fits(spotIndex).SpotFits3D = single(GaussFit);
    spotsFrame.Fits(spotIndex).SpotFits3DSE = single(FitDeltas);
    x1 = single(GaussFit(3) + xMin - 1);
    y1 = single(GaussFit(2) + yMin - 1);
    z1 = single(GaussFit(4) + zMin - 1);
    spotsFrame.Fits(spotIndex).GaussPos3D = single([x1,y1,z1]);
    spotsFrame.Fits(spotIndex).npFlag = isSymmetricPositiveDefinite;
    spotsFrame.Fits(spotIndex).int_dims = integrationDims;
    spotsFrame.Fits(spotIndex).GaussPos3DSE = single(FitDeltas(2:4));
    spotsFrame.Fits(spotIndex).gauss3DIntensity = single(GaussIntegral);
    spotsFrame.Fits(spotIndex).gauss3DIntensityRaw = single(GaussIntegralRaw);
    spotsFrame.Fits(spotIndex).gauss3DIntensitySE = single(GaussIntegralSE);
end
% [SpotsFr.Fits(spot).fits3D, SpotsFr.Fits(spot).gauss3DIntensity,...
%     SpotsFr.Fits(spot).fits3DCI95,...
%     SpotsFr.Fits(spot).gauss3DIntensityCI95] = fitGaussian3D(snip3D, initial_params, zStep, pixelSize, fitOptions{:});
% dx = SpotsFr.Fits(spot).fits3D(2) - SpotsFr.Fits(spot).fits3DCI95(2, 1);
% dy = SpotsFr.Fits(spot).fits3D(3) - SpotsFr.Fits(spot).fits3DCI95(3, 1);
% dz = SpotsFr.Fits(spot).fits3D(4) - SpotsFr.Fits(spot).fits3DCI95(4, 1);
% SpotsFr.Fits(spot).GaussPosCI95 = single([dx, dy, dz]);

%%%%%%%%%%%%%%%%
if displayFigures
    displayFigs(snip3D, GaussFit);
end

end

function displayFigs(snip3D, GaussFit)

close all;

[xMesh, yMesh, zMesh] = meshgrid(1:size(snip3D,1),1:size(snip3D,2),1:size(snip3D,3));
f1 = figure;
tiledlayout(f1,'flow')
nexttile;
imshow(imresize(max(snip3D,[],3), 10), [])
title('snippet max projection')
nexttile;
xDim = size(snip3D,1);
yDim = size(snip3D,2);
zDim = size(snip3D,3);
dimensionVector = [yDim, xDim, zDim];
gaussianFunc = @(params) simulate3DGaussianRho(dimensionVector, params) + params(11) + ...
    params(12)*xMesh + params(13)*yMesh + params(14)*zMesh;
gaussian = gaussianFunc(GaussFit);
isosurface(gaussian, 3);
axis tight

f2 = figure();
[yz,xz] = meshgrid(1:size(snip3D,2), 1:size(snip3D,1));
[yx,zx] = meshgrid(1:size(snip3D,2), 1:size(snip3D,3));
[xy,zy] = meshgrid(1:size(snip3D,1), 1:size(snip3D,3));
for z = 1:size(zx,1)
    maxxFit = max(gaussian, [], 1);
    maxyFit = max(gaussian, [], 2);
    maxxRaw = max(snip3D, [], 1);
    maxyRaw = max(snip3D, [], 2);
    projxFit(z,:) = maxxFit(:,:,z);
    projyFit(z,:) = maxyFit(:,:,z);
    projxRaw(z,:) = maxxRaw(:,:,z);
    projyRaw(z,:) = maxyRaw(:,:,z);
end
projz = max(gaussian, [], 3);
m = max(snip3D(:));
axesxyFit = axes(f2);
axesyzFit = axes(f2);
axesxzFit = axes(f2);
axesxyRaw = axes(f2);
axesyzRaw = axes(f2);
axesxzRaw = axes(f2);
surf(yz,xz, projz)
title(axesxyFit,{'Gaussian fit in X,Y';'(Z projected)'})
xlabel(axesxyFit,'x-axis')
ylabel(axesxyFit,'y-axis')

subplot(3,2,4,axesyzFit)
surf(yx,zx, projxFit)
title(axesyzFit,{'Gaussian fit in Y,Z';'(X projected)'})
xlabel(axesyzFit,'y-axis')
ylabel(axesyzFit,'z-axis')

subplot(3,2,6,axesxzFit)
surf(xy,zy, projyFit)
title(axesxzFit,{'Gaussian fit in X,Z';'(Y projected)'})
xlabel(axesxzFit,'x-axis')
ylabel(axesxzFit,'z-axis')

%%%% Raw Data Plotting %%%%
subplot(3,2,1,axesxyRaw)
surf(yz, xz, max(snip3D, [], 3))
title(axesxyRaw,{'Raw spot data in X,Y';'(Z projected)'})
xlabel(axesxyRaw,'x-axis')
ylabel(axesxyRaw,'y-axis')

subplot(3,2,3,axesyzRaw)
surf(yx, zx, projxRaw)
title(axesyzRaw,{'Raw spot data in Y,Z';'(X projected)'})
xlabel(axesyzRaw,'y-axis')
ylabel(axesyzRaw,'z-axis')

subplot(3,2,5,axesxzRaw)
surf(xy, zy, projyRaw)
title(axesxzRaw,{'Raw spot data in X,Z';'(Y projected)'})
xlabel(axesxzRaw,'x-axis')
ylabel(axesxzRaw,'z-axis')

zlim([0, m]);

% waitforbuttonpress;

end