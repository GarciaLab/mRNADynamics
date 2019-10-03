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
if nSpots == 2
    [GaussParams1, SpotInt1, GaussCI1, IntCI1, GaussParams2,SpotInt2,GaussCI2,IntCI2,offset, offsetCI] = ...
        fit3DGaussian2spot(snip3D,pixelSize);

    % spot 1 position
    SpotsFr.Fits(spot).Spot1Fits3D = single(GaussParams1);	
    x1 = single(GaussParams1(3) + min(xRange));
    y1 = single(GaussParams1(2) + min(yRange));
    z1 = single(GaussParams1(4) + min(zRange));
    SpotsFr.Fits(spot).Spot1Pos = single([x1,y1,z1]);
    SpotsFr.Fits(spot).Spot1GaussPosCI95 = single(GaussCI1(2:4,:));
    SpotsFr.Fits(spot).Spot1Int3D = single(SpotInt1);
    SpotsFr.Fits(spot).Spot1Fits3DCI95 = single(IntCI1);

    % spot 2 position
    SpotsFr.Fits(spot).Spot2Fits3D = single(GaussParams2);    
    x2 = single(GaussParams2(3) + min(xRange));
    y2 = single(GaussParams2(2) + min(yRange));
    z2 = single(GaussParams2(4) + min(zRange));
    SpotsFr.Fits(spot).Spot2Pos = single([x2,y2,z2]);
    SpotsFr.Fits(spot).Spot2GaussPosCI95 = single(GaussCI2(2:4,:));
    SpotsFr.Fits(spot).Spot2Int3D = single(SpotInt2);
    SpotsFr.Fits(spot).Spot2Fits3DCI95 = single(IntCI2);

    % combined metrics
    SpotsFr.Fits(spot).gauss3DIntensity = single(SpotInt1 + SpotInt2);
    SpotsFr.Fits(spot).gauss3DIntensityCI95 = single(sqrt(IntCI1^2+IntCI2^2));
    SpotsFr.Fits(spot).Offset3D = offset;
    SpotsFr.Fits(spot).OffsetCI = offsetCI;
    posArray = vertcat([x1,y1,z1],[x2,y2,z2]);
    SpotsFr.Fits(spot).GaussPos = (SpotsFr.Fits(spot).Spot1Int3D *posArray(1,:)...
                + SpotsFr.Fits(spot).Spot2Int3D *posArray(2,:)) /...
                (SpotsFr.Fits(spot).Spot1Int3D+SpotsFr.Fits(spot).Spot2Int3D);
    % calculate error in CM position
    PosDeltas1 = diff(GaussCI1(2:4,:)');
    PosDeltas2 = diff(GaussCI2(2:4,:)');
    PosDeltasTot = sqrt((SpotInt1*PosDeltas1).^2+(SpotInt2*PosDeltas2).^2);
    SpotsFr.Fits(spot).GaussPosCI95 = [SpotsFr.Fits(spot).GaussPos-PosDeltasTot 
        SpotsFr.Fits(spot).GaussPos+PosDeltasTot]';
    
elseif nSpots == 1      
    [GaussParams1, SpotIntegral, GaussCI, IntegralCI, ~] = fit3DGaussian(snip3D,pixelSize);
    SpotsFr.Fits(spot).SpotFits3D = single(GaussParams1);    
    x1 = single(GaussParams1(3) + min(xRange));
    y1 = single(GaussParams1(2) + min(yRange));
    z1 = single(GaussParams1(4) + min(zRange));
    SpotsFr.Fits(spot).GaussPos = single([x1,y1,z1]);
    SpotsFr.Fits(spot).GaussPosCI95 = single(GaussCI(2:4,:));
    SpotsFr.Fits(spot).gauss3DIntensity = single(SpotIntegral);
    SpotsFr.Fits(spot).gauss3DIntensityCI95 = single(IntegralCI);
end    
% [SpotsFr.Fits(spot).fits3D, SpotsFr.Fits(spot).gauss3DIntensity,...
%     SpotsFr.Fits(spot).fits3DCI95,...
%     SpotsFr.Fits(spot).gauss3DIntensityCI95] = fitGaussian3D(snip3D, initial_params, zStep, pixelSize, fitOptions{:});
% dx = SpotsFr.Fits(spot).fits3D(2) - SpotsFr.Fits(spot).fits3DCI95(2, 1);
% dy = SpotsFr.Fits(spot).fits3D(3) - SpotsFr.Fits(spot).fits3DCI95(3, 1);
% dz = SpotsFr.Fits(spot).fits3D(4) - SpotsFr.Fits(spot).fits3DCI95(4, 1);
% SpotsFr.Fits(spot).GaussPosCI95 = single([dx, dy, dz]);

end