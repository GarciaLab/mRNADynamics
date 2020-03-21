function SpotsFr = fitSnip3D(SpotsFr, spotChannel, spot, frame, Prefix,...
    PreProcPath, FrameInfo, nSpots)

%%


% extract basic fit parameters
s = SpotsFr.Fits(spot);
xSize = FrameInfo(1).PixelsPerLine;
ySize = FrameInfo(1).LinesPerFrame;
pixelSize = FrameInfo(1).PixelSize*1000; %nm
zStep = FrameInfo(1).ZStep*1000;
zMax = FrameInfo(1).NumberSlices+2;
snipDepth = uint8(ceil(2500/zStep));

% NL: Need to make this independent of 2D fit info 
bZ = s.brightestZ;
xSpot = s.xDoG(s.z==bZ);
ySpot = s.yDoG(s.z==bZ);

if isfield(s, 'snippet_size') && ~isempty(s.snippet_size)
    snippet_size = s.snippet_size;
else
    snippet_size = round(1500/pixelSize); % (in pixels)set to be around 1.5 um 
end
snippet_size = uint16(snippet_size(1));


zBot = max([1,bZ - snipDepth]);
zTop = min([zMax, bZ + snipDepth]);
zRange = zBot:zTop;
xRange = max([1,xSpot-snippet_size]):min([xSize,xSpot+snippet_size]);
yRange = max([1,ySpot-snippet_size]):min([ySize,ySpot+snippet_size]);
snip3D = NaN(numel(yRange),numel(xRange),numel(zBot:zTop));

%%
iter = 1;
for z = zRange    
    FullSlice=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(frame,3)...
        ,'_z' iIndex(z,2) '_ch' iIndex(spotChannel,2) '.tif']);
    snip3D(:,:,iter) = double(FullSlice(yRange,xRange)); 
    iter = iter + 1;
end

%%
xm = single(min(xRange));
ym = single(min(yRange));
zm = single(min(zRange));
% disable two-spot fits for now
if nSpots == 2
    warning('two spot fitting not currently supported. Performing single spot fit')
    nSpots = 1;
end
if nSpots == 2
    [GaussParams1, GaussParams2, offset, GaussIntVec, centroid_mean,...
        GaussSE1, GaussSE2, offsetSE, GaussIntSEVec, centroid_se] = ...
        fit3DGaussian2spot(snip3D,pixelSize);
    
    % spot 1 position
    SpotsFr.Fits(spot).Spot1Fits3D = single(GaussParams1);
    SpotsFr.Fits(spot).Spot1FitsSE3D = single(GaussSE1);
    x1 = single(GaussParams1(3) + xm - 1);
    y1 = single(GaussParams1(2) + ym - 1);
    z1 = single(GaussParams1(4) + zm - 1);
    SpotsFr.Fits(spot).Spot1Pos3D = single([x1,y1,z1]);
    SpotsFr.Fits(spot).Spot1Int3D = GaussIntVec(1);
    SpotsFr.Fits(spot).Spot1IntSE3D = GaussIntSEVec(1);

    % spot 2 position
    SpotsFr.Fits(spot).Spot2Fits3D = single(GaussParams2);
    SpotsFr.Fits(spot).Spot2CI3D = single(GaussSE2);
    x2 = single(GaussParams2(3) + xm - 1);
    y2 = single(GaussParams2(2) + ym - 1);
    z2 = single(GaussParams2(4) + zm - 1);
    SpotsFr.Fits(spot).Spot2Pos3D = single([x2,y2,z2]);
    SpotsFr.Fits(spot).Spot2Int3D = GaussIntVec(2);
    SpotsFr.Fits(spot).Spot2IntSE3D = GaussIntSEVec(2);

    % combined metrics
    SpotsFr.Fits(spot).gauss3DIntensity = GaussIntVec(3);
    SpotsFr.Fits(spot).gauss3DIntensitySE = GaussIntSEVec(3);
    SpotsFr.Fits(spot).Offset3D = offset;
    SpotsFr.Fits(spot).Offset3DSE = offsetSE;   
    SpotsFr.Fits(spot).GaussPos3D = centroid_mean;
    SpotsFr.Fits(spot).GaussPos3DSE = centroid_se;

% single spot fit    
elseif nSpots == 1      
    [GaussFit, FitDeltas, GaussIntegral, GaussIntegralSE, GaussIntegralRaw,np_flag,int_dims]  = ...
        fit3DGaussianRho(snip3D,[pixelSize zStep]);  
                                                
    SpotsFr.Fits(spot).SpotFits3D = single(GaussFit);  
    SpotsFr.Fits(spot).SpotFits3DSE = single(FitDeltas);  
    x1 = single(GaussFit(3) + xm - 1);
    y1 = single(GaussFit(2) + ym - 1);
    z1 = single(GaussFit(4) + zm - 1);
    SpotsFr.Fits(spot).GaussPos3D = single([x1,y1,z1]);
    SpotsFr.Fits(spot).npFlag = np_flag;
    SpotsFr.Fits(spot).int_dims = int_dims;
    SpotsFr.Fits(spot).GaussPos3DSE = single(FitDeltas(2:4));
    SpotsFr.Fits(spot).gauss3DIntensity = single(GaussIntegral);
    SpotsFr.Fits(spot).gauss3DIntensityRaw = single(GaussIntegralRaw);
    SpotsFr.Fits(spot).gauss3DIntensitySE = single(GaussIntegralSE);    
end    
% [SpotsFr.Fits(spot).fits3D, SpotsFr.Fits(spot).gauss3DIntensity,...
%     SpotsFr.Fits(spot).fits3DCI95,...
%     SpotsFr.Fits(spot).gauss3DIntensityCI95] = fitGaussian3D(snip3D, initial_params, zStep, pixelSize, fitOptions{:});
% dx = SpotsFr.Fits(spot).fits3D(2) - SpotsFr.Fits(spot).fits3DCI95(2, 1);
% dy = SpotsFr.Fits(spot).fits3D(3) - SpotsFr.Fits(spot).fits3DCI95(3, 1);
% dz = SpotsFr.Fits(spot).fits3D(4) - SpotsFr.Fits(spot).fits3DCI95(4, 1);
% SpotsFr.Fits(spot).GaussPosCI95 = single([dx, dy, dz]);

end