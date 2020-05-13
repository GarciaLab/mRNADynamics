function spotsFrame= fitSnip3D(spotsFrame, spotChannel,...
    spotIndex, frame, Prefix,~, ~, nSpots, imStack)

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
    imStack = movieMat(:, :, :, frame, spotChannel);
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
    
  snip3D = imStack(yRange, xRange);
  
else
    
    n = 1;
    for z = zRange    
        
        FullSlice=imread([preFolder, filesep,liveExperiment.Prefix,'_',iIndex(frame,3)...
            ,'_z' iIndex(z,2) '_ch' iIndex(spotChannel,2) '.tif']);
        
        snip3D(:,:,n) = double(FullSlice(yRange,xRange)); 
        n = n + 1;
        
    end
    
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

end