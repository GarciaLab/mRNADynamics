function SpotsFr = fitSnip3D(SpotsFr, spotChannel, spot, frame, Prefix, PreProcPath, ProcPath, FrameInfo, dogs, displayFigures, saveType)


s = SpotsFr.Fits(spot);
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
nameSuffix = ['_ch', iIndex(spotChannel, 2)];

dogProb = 'DOG_';

isZPadded = false;
try
firstdogname = [dogProb, Prefix, '_', iIndex(1, 3), '_z', iIndex(1, 2), nameSuffix];

firstdogpath = [ProcPath, filesep,Prefix,'_',filesep,'dogs',filesep,firstdogname];

matsPresent = exist([firstdogpath, '.mat'], 'file');
tifsPresent = exist([firstdogpath, '.tif'], 'file');

if ~strcmpi(saveType, 'none')
    if tifsPresent & ~matsPresent
        saveType = '.tif';
    elseif matsPresent & ~tifsPresent
        saveType = '.mat';
    elseif matsPresent & tifsPresent
        error('not sure which files to pick. check your processed folder.');
    end
end

firstdogpath = [firstdogpath, saveType];

if strcmpi(saveType, '.tif')
    try
        firstDoG = imread(firstdogpath);
    catch
        firstdogname = ['prob', Prefix, '_', iIndex(1, 3), '_z', iIndex(1, 2), nameSuffix];
        dogProb = 'prob';
        firstdogpath = [ProcPath, filesep,Prefix,'_',filesep,'dogs',filesep,firstdogname];
        firstdogpath = [firstdogpath, saveType];
        firstDoG = imread(firstdogpath);
        
    end
elseif strcmpi(saveType, '.mat')
    load(firstdogpath, 'plane');
    firstDoG = plane;
elseif strcmpi(saveType, 'none')
    firstDoG = dogs(:, :, 1, 1);
end

if sum(firstDoG(:)) == 0
    isZPadded = true;
end

end
FullDoGSlice = [];

%%
k = 1;
for z = zBot:zTop
    
    if isZPadded
        dogZ = z;
    else
        dogZ = z-1;
    end
    
    if z > 1 && z < zMax
        FullSlice=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(frame,3)...
            ,'_z' iIndex(z,2) '_ch' iIndex(spotChannel,2) '.tif']);
        try
            if isempty(dogs)
                dog_name = [dogProb, Prefix, '_', iIndex(frame, 3), '_z', iIndex(dogZ, 2), nameSuffix, saveType];
                dog_full_path = [ProcPath, filesep,Prefix,'_',filesep,'dogs',filesep,dog_name];
                if ~isempty(dir([ProcPath, filesep, Prefix,'_', filesep, 'dogs']))
                    if strcmpi(saveType, '.tif')
                        FullDoGSlice= double(imread(dog_full_path));
                    elseif strcmpi(saveType, '.mat')
                        dog_name = [dogProb, Prefix, '_', iIndex(frame, 3), '_z', iIndex(dogZ, 2), nameSuffix, saveType];
                        dog_full_path = [ProcPath, filesep,Prefix,'_',filesep,'dogs',filesep,dog_name];
                        load(dog_full_path, 'plane');
                        FullDoGSlice = double(plane);
                    end
                end
            else
                FullDoGSlice = dogs(:, :, dogZ, frame);
            end
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

[SpotsFr.Fits(spot).fits3D, SpotsFr.Fits(spot).gauss3DIntensity,...
    SpotsFr.Fits(spot).fits3DCI95,...
    SpotsFr.Fits(spot).gauss3DIntensityCI95] = fitGaussian3D(snip3D, initial_params, zStep, pixelSize, fitOptions{:});

%cast all as singles so the addition will work properly
snippet_size = single(snippet_size); xSpot = single(xSpot); ySpot = single(ySpot); snipDepth = single(snipDepth); bZ = single(bZ);

x = SpotsFr.Fits(spot).fits3D(2) - snippet_size + xSpot;
y = SpotsFr.Fits(spot).fits3D(3) - snippet_size + ySpot;
z = SpotsFr.Fits(spot).fits3D(4) - snipDepth + bZ;

% dxLow = SpotsFr.Fits(spot).fits3DCI95(2, 1) - snippet_size + xSpot;
% dxHigh = SpotsFr.Fits(spot).fits3DCI95(2, 2) - snippet_size + xSpot;
% dyLow = SpotsFr.Fits(spot).fits3DCI95(3, 1) - snippet_size + ySpot;
% dyHigh = SpotsFr.Fits(spot).fits3DCI95(3, 2) - snippet_size + ySpot;
% dzLow = SpotsFr.Fits(spot).fits3DCI95(4, 1) - snipDepth + bZ;
% dzHigh = SpotsFr.Fits(spot).fits3DCI95(4, 2) - snipDepth + bZ;

dx = x - SpotsFr.Fits(spot).fits3DCI95(2, 1);
dy = y - SpotsFr.Fits(spot).fits3DCI95(3, 1);
dz = z - SpotsFr.Fits(spot).fits3DCI95(4, 1);

SpotsFr.Fits(spot).GaussPos = single([x,y,z]);
% SpotsFr.Fits(spot).GaussPosCI95 = single([dxLow,dxHigh; dyLow, dyHigh; dzLow, dzHigh]);
SpotsFr.Fits(spot).GaussPosCI95 = single([dx, dy, dz]);

try
    if ~isempty(dogSnip3D)
        midind = ceil(size(dogSnip3D,3)/2);
        if size(dogSnip3D,3) > 2
            SpotsFr.Fits(spot).ampdog3 = single(sum(sum(sum(dogSnip3D(:,:,midind-1:midind+1)))) );
            SpotsFr.Fits(spot).ampdog3Max = single(max(max(max(dogSnip3D(:,:,midind-1:midind+1))) ));
        else
            SpotsFr.Fits(spot).ampdog3 = [];
            SpotsFr.Fits(spot).ampdog3Max = single(max(max(max(dogSnip3D(:,:,:)))));
        end
    end
end

%this is a flag that the fit was done over few z-frames so the
%user can decide if they want to keep the fit or not
if k < 3
    SpotsFr.Fits(spot).weeFit = true;
else
    SpotsFr.Fits(spot).weeFit = false;
end

end