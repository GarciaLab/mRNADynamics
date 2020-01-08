% Use metadata to stitc tiles together based on position info 


%Hardcoded variables for initial implementation 
clear all, close all
SourcePath = 'E:/Gabriella/LivemRNA\Data\RawDynamicsData';
Prefix = '2019-12-09-4xVasamEGFPVK22Homo-HbNbHomo-Anterior-27_5C-NC13NC14';
Date = '2019-12-09';
EmbryoName = '4xVasamEGFPVK22Homo-HbNbHomo-Anterior-27_5C-NC13NC14';
SurfTileDataFolder = [SourcePath, filesep, Date, filesep, EmbryoName,filesep,...
    'FullEmbryo\SurfFileStacks'];

%% 

LIFPath = [SourcePath, filesep, Date, filesep, EmbryoName,filesep,...
    'FullEmbryo\SurfTile.lif'];
NSlices = 4;
NChannels = 1;
fiducialChannel = 1;
framesIndex = 1;


r = bfGetReader();
% Decorate the reader with the Memoizer wrapper
r = loci.formats.Memoizer(r);
r.setId(LIFPath);
LIFImages = bfopen(LIFPath);
LIFMeta = LIFImages{:, 4};
r.close();

%HisSlices = generateHisSlices(images, NSlices, NChannels, fiducialChannel, framesIndex, seriesIndex)
TileIndex = 1;
T1ImageSlices = generateHisSlicesTileScan(LIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile1 = max(T1ImageSlices, [], 3);
filepath = [SourcePath, filesep, Date, filesep, EmbryoName,filesep,...
    'FullEmbryo\TestDir'];
imwrite(tile1, [filepath, filesep, 'tile1.tif']); 

TileIndex = 2;
T2ImageSlices = generateHisSlicesTileScan(LIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile2 = max(T2ImageSlices,[],3);
imwrite(tile2, [filepath, filesep, 'tile2.tif']); 

TileIndex = 3;
T3ImageSlices = generateHisSlicesTileScan(LIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile3 = max(T3ImageSlices, [], 3);
%imwrite(tile3, [filepath, filesep, 'tile3.tif']); 

TileIndex = 4;
T4ImageSlices = generateHisSlicesTileScan(LIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile4 = max(T4ImageSlices, [], 3);
%imwrite(tile4, [filepath, filesep, 'tile4.tif']);

%% 
PixelSize = .241; % units: um
PixelSize_m = PixelSize*10^(-6);
t1x = 6.02199091218657E-02; % units: m
t1y = 3.96804811614957E-02; % units: m
t2x = 6.02199091218657E-02; % units: m
t2y = 3.99019097329243E-02; % units: m
t3x = 6.04413376932943E-02; % units: m
t3y = 3.96804811614957E-02; % units: m
t4x = 6.04413376932943E-02;% units: m
t4y = 3.99019097329243E-02; % units: m


dy = (t2y-t1y)/PixelSize_m;
dx = (t3x-t1x)/PixelSize_m;

img = zeros(2048, 2048);
img(1:1024, 1:1024) = tile3;
img(31:1054, 920:1943) = tile4;
img(920:1943, 1:1024) = tile1;
img(950:2243, 920:1943) = tile2;

imshow(img)

%% 

MidLIFPath = [SourcePath, filesep, Date, filesep, EmbryoName,filesep,...
    'FullEmbryo\MidTile.lif'];
NSlices = 4;
NChannels = 1;
fiducialChannel = 1;
framesIndex = 1;


r = bfGetReader();
% Decorate the reader with the Memoizer wrapper
r = loci.formats.Memoizer(r);
r.setId(LIFPath);
MidLIFImages = bfopen(MidLIFPath);
LIFMeta = LIFImages{:, 4};
r.close();

%HisSlices = generateHisSlices(images, NSlices, NChannels, fiducialChannel, framesIndex, seriesIndex)
TileIndex = 1;
T1ImageSlices = generateHisSlicesTileScan(MidLIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile1 = max(T1ImageSlices, [], 3);
filepath = [SourcePath, filesep, Date, filesep, EmbryoName,filesep,...
    'FullEmbryo\MidTestDir'];
imwrite(tile1, [filepath, filesep, 'tile1.tif']); 

TileIndex = 2;
T2ImageSlices = generateHisSlicesTileScan(MidLIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile2 = max(T2ImageSlices,[],3);
imwrite(tile2, [filepath, filesep, 'tile2.tif']); 

TileIndex = 3;
T3ImageSlices = generateHisSlicesTileScan(MidLIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile3 = max(T3ImageSlices, [], 3);
imwrite(tile3, [filepath, filesep, 'tile3.tif']); 

TileIndex = 4;
T4ImageSlices = generateHisSlicesTileScan(MidLIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile4 = max(T4ImageSlices, [], 3);
imwrite(tile4, [filepath, filesep, 'tile4.tif']);



TileIndex = 5;
T5ImageSlices = generateHisSlicesTileScan(MidLIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile5 = max(T5ImageSlices, [], 3);
imwrite(tile5, [filepath, filesep, 'tile5.tif']);


TileIndex = 6;
T6ImageSlices = generateHisSlicesTileScan(MidLIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile6 = max(T6ImageSlices, [], 3);
imwrite(tile6, [filepath, filesep, 'tile6.tif']);


TileIndex = 7;
T7ImageSlices = generateHisSlicesTileScan(MidLIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile7 = max(T7ImageSlices, [], 3);
imwrite(tile7, [filepath, filesep, 'tile7.tif']);


TileIndex = 8;
T8ImageSlices = generateHisSlicesTileScan(MidLIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile8 = max(T8ImageSlices, [], 3);
imwrite(tile8, [filepath, filesep, 'tile8.tif']);


TileIndex = 9;
T9ImageSlices = generateHisSlicesTileScan(MidLIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile9 = max(T9ImageSlices, [], 3);
imwrite(tile9, [filepath, filesep, 'tile9.tif']);

% figure(1)
% imshow(tile1)
% figure(2)
% imshow(tile2)
% figure(3)
% imshow(tile3)
% figure(4)
% imshow(tile4)
% figure(5)
% imshow(tile5)
% figure(6)
% imshow(tile6)
% figure(7)
% imshow(tile7)
% figure(8)
% imshow(tile8)
% figure(9)
% imshow(tile9)

img = zeros(3072, 3072);
% img(1:1024,1:1024) = tile7;
% img(1:1024,1025:2048) = tile8;
% img(1:1024,2049:3072) = tile9;
% img(1025:2048,1:1024) = tile4;
% img(1025:2048,1025:2048) = tile5;
% img(1025:2048,2049:3072) = tile6;
% img(2049:3072, 1:1024) = tile1;
% img(2049:3072, 1025:2048) = tile2;
% img(2049:3072, 2049:3072) = tile3;
img(1:1024,1:1024) = tile7;
img(1:1024,920:1943) = tile8;
img(1:1024,1839:2862) = tile9;
img(920:1943,1:1024) = tile4;
img(920:1943,920:1943) = tile5;
img(920:1943,1839:2862) = tile6;
img(1839:2862, 1:1024) = tile1;
img(1839:2862, 920:1943) = tile2;
img(1839:2862, 1839:2862) = tile3;



imshow(img)


