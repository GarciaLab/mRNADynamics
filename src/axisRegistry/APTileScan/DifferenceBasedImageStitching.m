
%Hardcoded variables for initial implementation 
%clear all, close all
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
imwrite(tile3, [filepath, filesep, 'tile3.tif']); 

TileIndex = 4;
T4ImageSlices = generateHisSlicesTileScan(LIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile4 = max(T4ImageSlices, [], 3);
imwrite(tile4, [filepath, filesep, 'tile4.tif']);

%% Equivalent to normalizer function 
width = 80;
PixelSize = .241; % units: um
% This choice of sigma (= 10% of the nuclear diameter) keeps nuclei as
% features
sigma = .6/PixelSize;
fs = 5*sigma;
im1 = double(tile1);
im1gauss = imgaussfilt(im1, sigma);
im1g2 = im1-im1gauss;
im1g3 = (im1-im1gauss)./sqrt(im1gauss);

%im = im - imfilter(im, fspecial('average', width)); %mean normalize
%im = im ./ (sqrt(imfilter(im .* im, fspecial('average', width))));
%im = im / std(im(:));
%im2 = double(im); 
figure(1)
imshow(im1)
figure(2) 
imshow(im1gauss)
figure(3)
imshow(im1g2)
figure(4)
imshow(im1g3)

im2 = double(tile2);
im2gauss = imgaussfilt(im2, sigma);
im2g2 = im2-im2gauss;
im2g3 = (im2-im2gauss)./sqrt(im2gauss);

%im = im - imfilter(im, fspecial('average', width)); %mean normalize
%im = im ./ (sqrt(imfilter(im .* im, fspecial('average', width))));
%im = im / std(im(:));
%im2 = double(im); 
figure(5)
imshow(im2)
figure(6) 
imshow(im2gauss)
figure(7)
imshow(im2g2)
figure(8)
imshow(im2g3)
% 
% upsampled = sp.ndimage.zoom(img, 10, order=3)
%     conv_scaled = sp.ndimage.filters.gaussian_filter(upsampled, 10*sigma, truncate = 10*truncate)
%     coords_scaled = peak_local_max(conv_scaled, min_distance=min_distance*10, threshold_abs=0.0001)
%     if len(coords_scaled) > 0:
%         # THIS LINE IS PROBABLY AN ERROR! SHOULD BE sigma and truncate!!
%         #conv = sp.ndimage.filters.gaussian_filter(img, 10*sigma, truncate = 10*truncate)
%         conv = sp.ndimage.filters.gaussian_filter(img, sigma, truncate = truncate)
%         coords = coords_scaled/10.0
%         intensities = [conv[np.int(np.round(x[0])), np.int(np.round(x[1]))] for x in coords]
%         return coords, intensities
%     else:
%         return [], []



%% 
%Second set of images for sigma = 6 um
PixelSize = .241; % units: um
% This choice of sigma (= 10% of the nuclear diameter) keeps nuclei as
% features
sigma2 = 6/PixelSize;
fs = 5*sigma;
im1g4= imgaussfilt(im1, sigma2);
im1g5 = im1-im1g4;
im1g6 = (im1-im1g4)./sqrt(im1g4);
im1g7 = imbinarize(im1g4);


figure(1)
imshow(im1)
figure(2) 
imshow(im1g4)
figure(3)
imshow(im1g5)
figure(4)
imshow(im1g6)
figure(9)
imshow(im1g7)

im2g4 = imgaussfilt(im2, sigma2);
im2g5 = im2-im2g4;
im2g6 = (im2-im2g4)./sqrt(im2g4);
im2g7 = imbinarize(im2g4);



figure(5)
imshow(im2)
figure(6) 
imshow(im2g4)
figure(7)
imshow(im2g5)
figure(8)
imshow(im2g6)
figure(10)
imshow(im2g7)

im3 = double(tile3);
im3gauss = imgaussfilt(im3, sigma);
im3g2 = im3-im3gauss;
im3g3 = (im3-im3gauss)./sqrt(im3gauss);
im3g4 = imgaussfilt(im3, sigma2);
im3g5 = im3-im3g4;
im3g6 = (im3-im3g4)./sqrt(im3g4);
im3g7 = imbinarize(im3g4);

im4 = double(tile4);
im4gauss = imgaussfilt(im4, sigma);
im4g2 = im4-im4gauss;
im4g3 = (im4-im4gauss)./sqrt(im4gauss);
im4g4 = imgaussfilt(im4, sigma2);
im4g5 = im4-im4g4;
im4g6 = (im4-im4g4)./sqrt(im4g4);
im4g7 = imbinarize(im4g4);

x_min = 0; x_max = 110;
y_min = -110; y_max = 110;



[x0A, y0A, scores0A, areas0A] = autostitchLR(im1, im2, x_max, y_max, x_min, y_min); 
[x1A, y1A, scores1A, areas1A] = autostitchLR(im1gauss, im2gauss, x_max, y_max, x_min, y_min); 
[x2A, y2A, scores2A, areas2A] = autostitchLR(im1g2, im2g2, x_max, y_max, x_min, y_min); 
[x3A, y3A, scores3A, areas3A] = autostitchLR(im1g3, im2g3, x_max, y_max, x_min, y_min); 
[x4A, y4A, scores4A, areas4A] = autostitchLR(im1g4, im2g4, x_max, y_max, x_min, y_min); 
[x5A, y5A, scores5A, areas5A] = autostitchLR(im1g5, im2g5, x_max, y_max, x_min, y_min); 
[x6A, y6A, scores6A, areas6A] = autostitchLR(im1g6, im2g6, x_max, y_max, x_min, y_min); 
[x7A, y7A, scores7A, areas7A] = autostitchLR(im1g7, im2g7, x_max, y_max, x_min, y_min); 
%% 


[x0B, y0B, scores0B, areas0B] = autostitchLR(im3, im4, x_max, y_max, x_min, y_min); 
[x1B, y1B, scores1B, areas1B] = autostitchLR(im3gauss, im4gauss, x_max, y_max, x_min, y_min); 
[x2B, y2B, scores2B, areas2B] = autostitchLR(im3g2, im4g2, x_max, y_max, x_min, y_min); 
[x3B, y3B, scores3B, areas3B] = autostitchLR(im3g3, im4g3, x_max, y_max, x_min, y_min); 
[x4B, y4B, scores4B, areas4B] = autostitchLR(im3g4, im4g4, x_max, y_max, x_min, y_min); 
[x5B, y5B, scores5B, areas5B] = autostitchLR(im3g5, im4g5, x_max, y_max, x_min, y_min); 
[x6B, y6B, scores6B, areas6B] = autostitchLR(im3g6, im4g6, x_max, y_max, x_min, y_min); 
[x7B, y7B, scores7B, areas7B] = autostitchLR(im3g7, im4g7, x_max, y_max, x_min, y_min); 

% 
%% 

y_min = 0; x_min = -110;

[x0C, y0C, scores0C, areas0C] = autostitchBT(im1, im3, x_max, y_max, x_min, y_min); 
[x1C, y1C, scores1C, areas1C] = autostitchBT(im1gauss, im3gauss, x_max, y_max, x_min, y_min); 
[x2C, y2C, scores2C, areas2C] = autostitchBT(im1g2, im3g2, x_max, y_max, x_min, y_min); 
[x3C, y3C, scores3C, areas3C] = autostitchBT(im1g3, im3g3, x_max, y_max, x_min, y_min); 
[x4C, y4C, scores4C, areas4C] = autostitchBT(im1g4, im3g4, x_max, y_max, x_min, y_min); 
[x5C, y5C, scores5C, areas5C] = autostitchBT(im1g5, im3g5, x_max, y_max, x_min, y_min); 
[x6C, y6C, scores6C, areas6C] = autostitchBT(im1g6, im3g6, x_max, y_max, x_min, y_min); 
[x7C, y7C, scores7C, areas7C] = autostitchBT(im1g7, im3g7, x_max, y_max, x_min, y_min); 


%% 
y_min = 0; x_min = -110;

[x0D, y0D, scores0D, areas0D] = autostitchBT(im2, im4, x_max, y_max, x_min, y_min); 
[x1D, y1D, scores1D, areas1D] = autostitchBT(im2gauss, im4gauss, x_max, y_max, x_min, y_min); 
[x2D, y2D, scores2D, areas2D] = autostitchBT(im2g2, im4g2, x_max, y_max, x_min, y_min); 
[x3D, y3D, scores3D, areas3D] = autostitchBT(im2g3, im4g3, x_max, y_max, x_min, y_min); 
[x4D, y4D, scores4D, areas4D] = autostitchBT(im2g4, im4g4, x_max, y_max, x_min, y_min); 
[x5D, y5D, scores5D, areas5D] = autostitchBT(im2g5, im4g5, x_max, y_max, x_min, y_min); 
[x6D, y6D, scores6D, areas6D] = autostitchBT(im2g6, im4g6, x_max, y_max, x_min, y_min); 
[x7D, y7D, scores7D, areas7D] = autostitchBT(im2g7, im4g7, x_max, y_max, x_min, y_min); 



%% 

figure(1)
plot(areas0A, scores0A, '.')
title('scores0A')
figure(2)
plot(areas1A, scores1A, '.')
figure(3)
plot(areas2A, scores2A, '.')
figure(4)
plot(areas3A, scores3A, '.')
figure(5)
plot(areas4A, score42A, '.')
figure(6)
plot(areas5A, scores5A, '.')
figure(7)
plot(areas6A, scores6A, '.')
figure(8)
plot(areas7A, scores7A, '.')
% 
% fig, ax = plt.subplots()
% p = ax.imshow(f.T[::-1], cmap='viridis')
% 
% ax.set_xticks(np.arange(len(MLabels))[::2])
% ax.set_xticklabels(MLabels[::2])
% ax.set_yticks(np.arange(len(rLabels))[::2])
% ax.set_yticklabels(rLabels[::-2])
% ax.set_ylabel('Fluorescence Ratio')
% ax.set_xlabel('Maximum $[Hb_{total}]$ (nM)')
% 
% #plt.imshow(f, cmap = 'viridis')
% fig.colorbar(p)
% 
% plt.show()
%% 


% %% 
% [h, w] = size(im1);
% 
% 
% y_min = -150; y_max = 150;
% x_min = 0; x_max = 150;
% 
% XRange = x_min:1:x_max;
% YRange = y_min:1:y_max;
% scores0 = zeros(length(XRange), length(YRange));
% scores1 = zeros(length(XRange), length(YRange));
% scores2 = zeros(length(XRange), length(YRange));
% scores3 = zeros(length(XRange), length(YRange));
% scores4 = zeros(length(XRange), length(YRange));
% scores5 = zeros(length(XRange), length(YRange));
% scores6 = zeros(length(XRange), length(YRange));
% scores7 = zeros(length(XRange), length(YRange));
% for i=1:length(XRange)
%     x = XRange(i);
%     for j=1:length(YRange)
%         y = YRange(j);
%         [h, w] = size(im1);
%         if y>=0
%             lw0 = im1(y+1:h,w-x+1:w);
%             rw0 = im2(1:h-y, 1:x);
%             lw1 = im1gauss(y+1:h,w-x+1:w);
%             rw1 = im2gauss(1:h-y, 1:x);
%             lw2 = im1g2(y+1:h,w-x+1:w);
%             rw2 = im2g2(1:h-y, 1:x);
%             lw3 = im1g3(y+1:h,w-x+1:w);
%             rw3 = im2g3(1:h-y, 1:x);
%             lw4 = im1g4(y+1:h,w-x+1:w);
%             rw4 = im2g4(1:h-y, 1:x);
%             lw5 = im1g5(y+1:h,w-x+1:w);
%             rw5 = im2g5(1:h-y, 1:x);
%             lw6 = im1g6(y+1:h,w-x+1:w);
%             rw6 = im2g6(1:h-y, 1:x);
%             lw7 = im1g7(y+1:h,w-x+1:w);
%             rw7 = im2g7(1:h-y, 1:x);
%             area = (h-y)*x;
%         else
%             lw0 = im1(1:h+y,w-x+1:w);
%             rw0 = im2(-y+1:h,1:x);
%             lw1 = im1gauss(1:h+y,w-x+1:w);
%             rw1 = im2gauss(-y+1:h,1:x);
%             lw2 = im1g2(1:h+y,w-x+1:w);
%             rw2 = im2g2(-y+1:h,1:x);
%             lw3 = im1g3(1:h+y,w-x+1:w);
%             rw3 = im2g3(-y+1:h,1:x);
%             lw4 = im1g4(1:h+y,w-x+1:w);
%             rw4 = im2g4(-y+1:h,1:x);
%             lw5 = im1g5(1:h+y,w-x+1:w);
%             rw5 = im2g5(-y+1:h,1:x);
%             lw6 = im1g6(1:h+y,w-x+1:w);
%             rw6 = im2g6(-y+1:h,1:x);
%             lw7 = im1g7(1:h+y,w-x+1:w);
%             rw7 = im2g7(-y+1:h,1:x);
%             area = (h+y)*x;
%         end
%         diffSquared0 = (lw0 - rw0).^ 2;
%         scores0(i,j) = sqrt(sum(diffSquared0(:)))/area;
%         diffSquared1 = (lw1 - rw1).^ 2;
%         scores1(i,j) = sqrt(sum(diffSquared1(:)))/area;
%         diffSquared2 = (lw2 - rw2).^ 2;
%         scores2(i,j) = sqrt(sum(diffSquared2(:)))/area;
%         diffSquared3 = (lw3 - rw3).^ 2;
%         scores3(i,j) = sqrt(sum(diffSquared3(:)))/area;
%         diffSquared4 = (lw4 - rw4).^ 2;
%         scores4(i,j) = sqrt(sum(diffSquared4(:)))/area;
%         diffSquared5 = (lw5 - rw5) .^ 2;
%         scores5(i,j) = sqrt(sum(diffSquared5(:)))/area;
%         diffSquared6 = (lw6 - rw6) .^ 2;
%         scores6(i,j) = sqrt(sum(diffSquared6(:)))/area;
%         diffSquared7 = (lw7 - rw7) .^ 2;
%         scores7(i,j) = sqrt(sum(diffSquared7(:)))/area;
%     end
% end
%% 

% 
% [x0ind,y0ind]=find(scores0==nanmin(nanmin(scores0)));
% x0 = XRange(x0ind);
% y0 = YRange(y0ind);
% 
% [x1ind,y1ind]=find(scores1==nanmin(nanmin(scores1)));
% x1 = XRange(x1ind);
% y1 = YRange(y1ind);
% 
% [x2ind,y2ind]=find(scores2==nanmin(nanmin(scores2)));
% x2 = XRange(x2ind);
% y2 = YRange(y2ind);
% 
% [x3ind,y3ind]=find(scores3==nanmin(nanmin(scores3)));
% x3 = XRange(x3ind);
% y3 = YRange(y3ind);
% 
% 
% [x4ind,y4ind]=find(scores4==nanmin(nanmin(scores4)));
% x4 = XRange(x4ind);
% y4 = YRange(y4ind);
% 
% 
% [x5ind,y5ind]=find(scores4==nanmin(nanmin(scores5)));
% x5 = XRange(x5ind);
% y5 = YRange(y5ind);
% 
% 
% [x6ind,y6ind]=find(scores6==nanmin(nanmin(scores6)));
% x6 = XRange(x6ind);
% y6 = YRange(y6ind);
% 
% 
% [x7ind,y7ind]=find(scores7==nanmin(nanmin(scores7)));
% x7 = XRange(x7ind);
% y7 = YRange(y7ind);
% 
% 
% %% Tiles 3 and 4
% PixelSize = .241; % units: um
% % This choice of sigma (= 10% of the nuclear diameter) keeps nuclei as
% % features
% sigma = .6/PixelSize;
% fs = 5*sigma;
% im3 = double(tile3);
% im3gauss = imgaussfilt(im3, sigma);
% im3g2 = im3-im3gauss;
% im3g3 = (im3-im3gauss)./sqrt(im3gauss);
% 
% %im = im - imfilter(im, fspecial('average', width)); %mean normalize
% %im = im ./ (sqrt(imfilter(im .* im, fspecial('average', width))));
% %im = im / std(im(:));
% %im2 = double(im); 
% figure(1)
% imshow(im3)
% figure(2) 
% imshow(im3gauss)
% figure(3)
% imshow(im3g2)
% figure(4)
% imshow(im3g3)
% 
% im2 = double(tile2);
% im2gauss = imgaussfilt(im2, sigma);
% im2g2 = im2-im2gauss;
% im2g3 = (im2-im2gauss)./sqrt(im2gauss);
% 
% %im = im - imfilter(im, fspecial('average', width)); %mean normalize
% %im = im ./ (sqrt(imfilter(im .* im, fspecial('average', width))));
% %im = im / std(im(:));
% %im2 = double(im); 
% figure(5)
% imshow(im2)
% figure(6) 
% imshow(im2gauss)
% figure(7)
% imshow(im2g2)
% figure(8)
% imshow(im2g3)
% % 
% % upsampled = sp.ndimage.zoom(img, 10, order=3)
% %     conv_scaled = sp.ndimage.filters.gaussian_filter(upsampled, 10*sigma, truncate = 10*truncate)
% %     coords_scaled = peak_local_max(conv_scaled, min_distance=min_distance*10, threshold_abs=0.0001)
% %     if len(coords_scaled) > 0:
% %         # THIS LINE IS PROBABLY AN ERROR! SHOULD BE sigma and truncate!!
% %         #conv = sp.ndimage.filters.gaussian_filter(img, 10*sigma, truncate = 10*truncate)
% %         conv = sp.ndimage.filters.gaussian_filter(img, sigma, truncate = truncate)
% %         coords = coords_scaled/10.0
% %         intensities = [conv[np.int(np.round(x[0])), np.int(np.round(x[1]))] for x in coords]
% %         return coords, intensities
% %     else:
% %         return [], []
% 
% 
% 
% %% 
% %Second set of images for sigma = 6 um
% PixelSize = .241; % units: um
% % This choice of sigma (= 10% of the nuclear diameter) keeps nuclei as
% % features
% sigma2 = 6/PixelSize;
% fs = 5*sigma;
% im3g4= imgaussfilt(im3, sigma2);
% im3g5 = im3-im3g4;
% im3g6 = (im3-im3g4)./sqrt(im3g4);
% im3g7 = imbinarize(im3g4);
% 
% 
% figure(1)
% imshow(im3)
% figure(2) 
% imshow(im3g4)
% figure(3)
% imshow(im3g5)
% figure(4)
% imshow(im3g6)
% figure(9)
% imshow(im3g7)
% 
% im2g4 = imgaussfilt(im2, sigma2);
% im2g5 = im2-im2g4;
% im2g6 = (im2-im2g4)./sqrt(im2g4);
% im2g7 = imbinarize(im2g4);
% 
% 
% 
% figure(5)
% imshow(im2)
% figure(6) 
% imshow(im2g4)
% figure(7)
% imshow(im2g5)
% figure(8)
% imshow(im2g6)
% figure(10)
% imshow(im2g7)


%% 
% figure(1)
% imshow(imstitchLR(im1, im2, x0, y0, [1 2]))
% figure(2)
% imshow(imstitchLR(im1, im2, x1, y1, [1 2]))
% figure(3)
% imshow(imstitchLR(im1, im2, x2, y2, [1 2]))
% %figure(4)
% %imshow(imstitchLR(im1, im2, x3, y3, [1 2]))
% figure(5)
% imshow(imstitchLR(im1, im2, x4, y4, [1 2]))
% %figure(6)
% %imshow(imstitchLR(im1, im2, x5, y5, [1 2]))
% figure(7)
% imshow(imstitchLR(im1, im2, x6, y6, [1 2]))
% figure(8)
% imshow(imstitchLR(im1, im2, x7, y7, [1 2]))


%% 



function [r_offset, c_offset, scores, areas] = autostitchLR(left, right, ro_max, co_max, ro_min, co_min) 
    % goes through where the widths are close, and match the two normalized
    % images to see how good they match
    if co_min < 0
        error(['column offset cannot be negative: tiled images must have an overlap',...
        'greater than or equal to zero. Choose a non-negative value for co_min.' ]);
    end
    [h, w] = size(left);
    RoRange = ro_min:1:ro_max;
    CoRange = co_min:1:co_max;
    scores = zeros(length(RoRange), length(CoRange));
    areas = zeros(length(RoRange), length(CoRange));
    for i=1:length(RoRange)
        ro = RoRange(i);
        for j=1:length(CoRange)
            co = CoRange(j);
            if ro >=0
                leftWindow = left(ro+1:h,w-co+1:co);
                rightWindow = right(1:h-ro, 1:co); 
                areas(i,j) = (h-ro)*co;
            else
                leftWindow = left(1:h+ro,w-co+1:w);
                rightWindow = right(-ro+1:h,1:co);
                areas(i,j) = (h+ro)*co;
            end

            diffSquared = (leftWindow - rightWindow) .^ 2;
            scores(i,j) = sqrt(sum(diffSquared(:)))/areas(i,j);
                
        end
    end
    
    [ro_index, co_index] = find(scores==nanmin(nanmin((scores))));
    
    r_offset = RoRange(ro_index);
    c_offset = YRange(co_index);
end

function [r_offset, c_offset, scores, areas] = autostitchBT(bottom, top, ro_max, co_max, ro_min, co_min) 
    % goes through where the widths are close, and match the two normalized
    % images to see how good they match
    if ro_min < 0
        error(['row offset cannot be negative: tiled images must have an overlap',...
        'greater than or equal to zero. Choose a non-negative value for ro_min.' ]);
    end
    [h, w] = size(bottom);
    RoRange = ro_min:1:ro_max;
    CoRange = co_min:1:co_max;
    scores = zeros(length(RoRange), length(CoRange));
    areas = zeros(length(RoRange), length(CoRange));
    for i=1:length(RoRange)
        ro = RoRange(i);
        for j=1:length(CoRange)
            co = CoRange(j);
            if co >=0
                topWindow = top(h-ro+1:h,co+1:w);
                bottomWindow = bottom(1:ro,1:w-co);
                areas(i,j) = (w-co)*ro;
            else
                topWindow = top(h-ro+1:h,1:w+co);
                bottomWindow = bottom(1:ro, -co+1:w);
                areas(i,j) = (w+co)*ro;
            end

            diffSquared = (topWindow - bottomWindow) .^ 2;
            scores(i,j) = sqrt(sum(diffSquared(:)))/areas(i,j);
                
        end
    end
    
    [ro_index, co_index] = find(scores==nanmin(nanmin((scores))));
    
    r_offset = RoRange(ro_index);
    c_offset = CoRange(co_index);
end





function imm2 = imstitchLR(left,right, xo1, yo1,order)
    %% stitches images given offsets. you specify the drawing order
    % [xo1 yo1] = offsets
    % order = draw order {1=left, 2 = right}. drawn left to right
    
    [h, w] = size(left);
    imm2 = zeros(h+abs(yo1), 2*w);

    % ?????
        for i=order
            switch i
                case 1
                    if yo1 >= 0
                      imm2(1:h, 1:w) = left;
                    else
                      imm2((abs(yo1)+1):(h+abs(yo1)), 1:w)=left;
                    end
                case 2
                    if yo1 >= 0
                        imm2((yo1 + 1):(yo1 + h), w+1-xo1:2*w-xo1) = right;
                    else
                        imm2(1:h, w+1-xo1:2*w-xo1) = right;
                    end
            end
        end
end

function imm2 = imstitchBT(bottom,top, xo1, yo1,order)
    %% stitches images given offsets. you specify the drawing order
    % [xo1 yo1] = offsets
    % order = draw order {1=top, 2 = bottom}. drawn top to bottom
    
    [h, w] = size(top);
    imm2 = zeros(2*h, w+abs(xo1));

    for i=order
        switch i
            case 1
                if xo1 >= 0
                  imm2(1:h, 1:w) = top;
                else
                  imm2(1:h, (abs(xo1)+1):(abs(xo1) + w)) = top;
                end
            case 2
                if xo1 >= 0
                    imm2(h+1-yo1:2*h-yo1, xo1+1:xo1+w) = bottom;
                    %imm2((yo1 + 1):(yo1 + size(left, 1)), w+1-xo1:2*w-xo1) = right;
                else
                    imm2(h+1-yo1:2*h-yo1, 1:w) = bottom;
                end
        end
    end
end




