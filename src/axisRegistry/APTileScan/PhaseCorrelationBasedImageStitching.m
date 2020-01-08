% Doesn't work 
% Based on the following link:
% https://thilinasameera.wordpress.com/2012/03/24/translation-invariant-image-registration-using-phase-correlation-panorama-imaging-on-matlab/


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
PixelSize = 241; % units: nm

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

%% Based on website 
% Define pararameters for code
image_width = 1024;
image_height = 1024;
image_left = tile1;
image_right = tile2;
%% 


%
hanning1D_row = hanning(image_width);
hanning1D_col = hanning(image_height);
hanning2D = hanning1D_row*hanning1D_col';

image_left_flt = hanning2D.*image_left;
image_right_flt = hanning2D.*image_right;

%% 

FFT_L = fft2(image_left_flt );
FFT_R = fft2(image_right_flt );
IMF = FFT_L.*conj(FFT_R);
CPS = IMF./abs(IMF);
%% 
im1 = image_left;
im2 = image_right;

F1 = fftshift(fft2(im1));
F2 = fftshift(fft2(im2));
 
[r1 c1] = size(im1);
[r2 c2] = size(im2);
 
IMF = F1.*conj(F2);
CPS = IMF./abs(IMF);
mag = (ifft2(F));
 
[X Y] = find((mag == (max(max(mag)))));
 
if im1(1,1)==im2(end-Y,end-X)
    img = zeros(Y+r2,X+c2);
    %img=cast(img,'uint8');
    img(1:r1,1:c1)=im1;
    img(Y:Y+r2-1,X:X+c2-1)=im2;
else
    img = zeros(Y+r2,X+c2);
    %img=cast(img,'uint8');
    img(1:r1,1:c1)=im1;
    img(Y:Y+r2-1,X:X+c2-1)=im2;
end
 

imshow(img);
%imout = imcrop(img);
