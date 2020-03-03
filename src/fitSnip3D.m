function SpotsFr = fitSnip3D(SpotsFr, spotChannel, spot, frame, Prefix,...
    PreProcPath, FrameInfo, nSpots, movieMat)

%%
% check for pre-existing 3D fit info
% prev_flag = 0;

% if ~isempty(SpotsFr(i).Fits)
%     fnames = fieldnames(SpotsFr(i).Fits);
%     prev_indices = contains(fnames,'3D');
%     SpotsFr(i).Fits = rmfield(SpotsFr(i).Fits,fnames(prev_indices));
%     if any(prev_indices) && prev_flag == 1
%         warning('previous 3D fit info detected. Removing...')
%     end
% end


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


snip3D = double(squeeze(movieMat(yRange,xRange, zRange, frame, spotChannel))); 


%%
xm = single(min(xRange));
ym = single(min(yRange));
zm = single(min(zRange));
if nSpots == 2
    [GaussParams1, GaussParams2, offset, GaussIntVec,...
        centroid_mean, GaussSE1, GaussSE2, offsetSE, GaussIntSEVec, centroid_se] = ...
    ...    
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
            
elseif nSpots == 1      
    [GaussFit, FitDeltas, GaussIntegral, GaussIntegralSE, GaussIntegralRaw]  = ...
        fit3DGaussian(snip3D,[pixelSize zStep]);
    SpotsFr.Fits(spot).SpotFits3D = single(GaussFit);  
    SpotsFr.Fits(spot).SpotFits3DSE = single(FitDeltas);  
    x1 = single(GaussFit(3) + xm - 1);
    y1 = single(GaussFit(2) + ym - 1);
    z1 = single(GaussFit(4) + zm - 1);
    SpotsFr.Fits(spot).GaussPos3D = single([x1,y1,z1]);
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