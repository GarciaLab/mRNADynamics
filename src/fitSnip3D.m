function SpotsCh = fitSnip3D(SpotsCh, channel, spot, frame, Prefix, PreProcPath, ProcPath, FrameInfo, dogs, displayFigures, saveType)

s = SpotsCh(frame).Fits(spot);
xSize = FrameInfo(1).PixelsPerLine;
ySize = FrameInfo(1).LinesPerFrame;
pixelSize = FrameInfo(1).PixelSize*1000; %nm
zStep = FrameInfo(1).ZStep;
zMax = FrameInfo(1).NumberSlices+2;

zoomFactor = 1; %replace this with zStep from FrameInfo later- AR 1/31/2019
snipDepth = uint8(round(3*zoomFactor));


bZ = s.brightestZ;
xSpot = s.xDoG(s.z==bZ);
ySpot = s.yDoG(s.z==bZ);

if isfield(s, 'snippet_size') && ~isempty(s.snippet_size)
    snippet_size = s.snippet_size;
else
    snippet_size = 13; %pixels
end
snippet_size = uint16(snippet_size(1));


zBot = bZ - snipDepth;
zTop = bZ + snipDepth;
snip3D = [];
dogSnip3D = [];
nameSuffix = ['_ch', iIndex(channel, 2)];


isZPadded = false;
  firstdogname = ['DOG_', Prefix, '_', iIndex(1, 3), '_z', iIndex(1, 2), nameSuffix, saveType];
     firstdogpath = [ProcPath, filesep,Prefix,'_',filesep,'dogs',filesep,firstdogname];

if strcmpi(saveType, '.tif')
   firstDoG = imread(firstdogpath);
elseif strcmpi(saveType, '.mat')
     load(firstdogpath, 'plane');
    firstDoG = plane;
elseif strcmpi(saveType, 'none')
    firstDoG = dogs(:, :, 1, 1);
end

if sum(firstDoG(:)) == 0
    isZPadded = true;
end


%%
k = 1;
for z = zBot:zTop
    
    if isZPadded
        dogZ = z;
    else
        dogZ = z-1;
        dogZMax = zMax - 2;
    end
    
    if z > 1 && z < zMax
        FullSlice=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(frame,3)...
            ,'_z' iIndex(z,2) '_ch' iIndex(channel,2) '.tif']);
        if isempty(dogs)
            dog_name = ['DOG_', Prefix, '_', iIndex(frame, 3), '_z', iIndex(dogZ, 2), nameSuffix, saveType];
            dog_full_path = [ProcPath, filesep,Prefix,'_',filesep,'dogs',filesep,dog_name];
            if ~isempty(dir([ProcPath, filesep, Prefix,'_', filesep, 'dogs']))
                if strcmpi(saveType, '.tif')
                    FullDoGSlice= double(imread(dog_full_path));
                elseif strcmpi(saveType, '.mat')
                    dog_name = ['DOG_', Prefix, '_', iIndex(frame, 3), '_z', iIndex(dogZ, 2), nameSuffix, saveType];
                    dog_full_path = [ProcPath, filesep,Prefix,'_',filesep,'dogs',filesep,dog_name];
                    load(dog_full_path, 'plane');
                    FullDoGSlice = double(plane);
                end
            else
                FullDoGSlice = [];
            end
        else
            FullDoGSlice = dogs(:, :, dogZ, frame);
        end
        snip3D(:,:,k) = double(FullSlice(max(1,ySpot-snippet_size):min(ySize,ySpot+snippet_size),...
            max(1,xSpot-snippet_size):min(xSize,xSpot+snippet_size))); %#ok<*SAGROW>
        if ~isempty(FullDoGSlice)
            dogSnip3D(:, :, k) = double(FullDoGSlice(max(1,ySpot-snippet_size):min(ySize,ySpot+snippet_size),...
                max(1,xSpot-snippet_size):min(xSize,xSpot+snippet_size))); %#ok<*SAGROW>
        else
            dogSnip3D = [];
        end
        k = k + 1;
    end
end
%%

%%
initial_params = [max(snip3D(:)), snipDepth + 1, 200/pixelSize , nanmean(s.Offset)]; %nm. empirically determined and seems to work width of spot psf
fitOptions = {};
if displayFigures
    fitOptions = [fitOptions,'displayFigures'];
end

[SpotsCh(frame).Fits(spot).fits3D, SpotsCh(frame).Fits(spot).gauss3DIntensity,...
    SpotsCh(frame).Fits(spot).fits3DCI95, SpotsCh(frame).Fits(spot).gauss3DIntensityCI95] = fitGaussian3D(snip3D, initial_params, zStep, pixelSize, fitOptions{:});

%cast all as singles so the addition will work properly
snippet_size = single(snippet_size); xSpot = single(xSpot); ySpot = single(ySpot); snipDepth = single(snipDepth); bZ = single(bZ);

x = SpotsCh(frame).Fits(spot).fits3D(2)  -snippet_size + xSpot;
y = SpotsCh(frame).Fits(spot).fits3D(3)  -snippet_size + ySpot;
z = SpotsCh(frame).Fits(spot).fits3D(4) -  snipDepth + bZ;

dxLow = SpotsCh(frame).Fits(spot).fits3DCI95(2, 1) - snippet_size + xSpot;
dxHigh = SpotsCh(frame).Fits(spot).fits3DCI95(2, 2) - snippet_size + xSpot;
dyLow = SpotsCh(frame).Fits(spot).fits3DCI95(3, 1) - snippet_size + ySpot;
dyHigh = SpotsCh(frame).Fits(spot).fits3DCI95(3, 2) - snippet_size + ySpot;
dzLow = SpotsCh(frame).Fits(spot).fits3DCI95(4, 1) - snipDepth + bZ;
dzHigh = SpotsCh(frame).Fits(spot).fits3DCI95(4, 2) - snipDepth + bZ;

SpotsCh(frame).Fits(spot).GaussPos = single([x,y,z]);
SpotsCh(frame).Fits(spot).GaussPosCI95 = single([dxLow,dxHigh; dyLow, dyHigh; dzLow, dzHigh]);

if ~isempty(dogSnip3D)
    midind = ceil(size(dogSnip3D,3)/2);
    if size(dogSnip3D,3) > 2
        SpotsCh(frame).Fits(spot).ampdog3 = single(sum(sum(sum(dogSnip3D(:,:,midind-1:midind+1)))) );
        SpotsCh(frame).Fits(spot).ampdog3Max = single(max(max(max(dogSnip3D(:,:,midind-1:midind+1))) ));
    else
        SpotsCh(frame).Fits(spot).ampdog3 = [];
        SpotsCh(frame).Fits(spot).ampdog3Max = single(max(max(max(dogSnip3D(:,:,:)))));
    end
end

%this is a flag that the fit was done over few z-frames so the
%user can decide if they want to keep the fit or not
if k < 3
    SpotsCh(frame).Fits(spot).weeFit = true;
else
    SpotsCh(frame).Fits(spot).weeFit = false;
end

end