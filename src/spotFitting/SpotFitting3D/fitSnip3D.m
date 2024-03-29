function [spotsFrame, fitInfo] = fitSnip3D(spotsFrame, spotIndex, liveExperiment, imStack, spotDims, nSpots)

FrameInfo = getFrameInfo(liveExperiment);

% extract basic fit parameters
spot = spotsFrame.Fits(spotIndex);
xSize = FrameInfo(1).PixelsPerLine;
ySize = FrameInfo(1).LinesPerFrame;
pixelSize_nm = FrameInfo(1).PixelSize*1000; %nm
zStep_nm = FrameInfo(1).ZStep*1000; %nm
zMax = FrameInfo(1).NumberSlices+1;
snipDepth = ceil(3500/zStep_nm); % sets Z depth above and below brightest plane that is extracted for fitting

% NL: Ideally we would make this independent of 2D fit info 
brightestZPlane = spot.brightestZ;
xSpot = spot.xDoG(spot.z==brightestZPlane);
ySpot = spot.yDoG(spot.z==brightestZPlane);

if isfield(spot, 'snippet_size') && ~isempty(spot.snippet_size)
    snippet_size = spot.snippet_size;
else
    snippet_size = round(1500/pixelSize_nm); % (in pixels)set to be around 1.5 um 
end
snippet_size = uint16(snippet_size(1));

%% %%%%%%%%%%%%%%%%%%%%%% generate 3D snip %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zBot = max([2,brightestZPlane - snipDepth]);
zTop = min([zMax-1, brightestZPlane + snipDepth]);
zRange = zBot:zTop;
xRange = max([1,xSpot-snippet_size]):min([xSize,xSpot+snippet_size]);
yRange = max([1,ySpot-snippet_size]):min([ySize,ySpot+snippet_size]);

snip3D = imStack(yRange, xRange, zRange);
  
%% %%%%%%%%%%%%%%%%%%%%%% perform fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get stack origins
xMin = single(min(xRange));
yMin = single(min(yRange));
zMin = single(min(zRange));

% Perform fits
fitInfo = fit3DGaussians(snip3D, pixelSize_nm, zStep_nm, spotDims, nSpots, []);

% add fit info
spotsFrame.Fits(spotIndex).SpotFitInfo3D = fitInfo;%single(GaussParams1);

% combined metrics
spotsFrame.Fits(spotIndex).GaussIntensity3D = fitInfo.GaussIntegralTot;
spotsFrame.Fits(spotIndex).GaussIntensity3DSE = fitInfo.GaussIntegralSETot;
spotsFrame.Fits(spotIndex).Offset3D = fitInfo.offset;
spotsFrame.Fits(spotIndex).GaussIntensity3DRaw = fitInfo.GaussIntegralRaw;
spotsFrame.Fits(spotIndex).GaussIntensity3DRawSE = fitInfo.GaussIntegralRawSE;

% spotsFrame.Fits(spotIndex).Offset3DSE = offsetSE;   
spotsFrame.Fits(spotIndex).GaussPos3D = fitInfo.SpotCentroid + [yMin xMin zMin] - 1.0;
spotsFrame.Fits(spotIndex).GaussPos3DSE = fitInfo.SpotCentroidSE;


end