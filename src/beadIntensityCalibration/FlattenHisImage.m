clear all, close all
Prefix='2019-12-10-4xVasamEGFPVK22Homo-HbNbHomo-Anterior-27_5C-NC14Peak';
dataFolder = ['E:\Gabriella\LivemRNA\Data\PreProcessedData', filesep, Prefix];
rootFolder = 'E:\Gabriella\LivemRNA\Data\RawDynamicsData\2019-12-10\4xVasamEGFPVK22Homo-HbNbHomo-Anterior-27_5C-NC14Peak';
D = dir([dataFolder, filesep, '*His*.tif']);
names = cell(length(D), 1);
for i=1:length(D)
    names{i}=[dataFolder,filesep,D(i).name];
end

rawDataDir = dir([rootFolder, filesep, '*2019-12-10*.lif']);
im2 = imread([D(1).folder, filesep, D(1).name]);
LIFPath = [rawDataDir(1).folder,filesep, rawDataDir(1).name];
DropboxFolder = 'E:\Gabriella\LivemRNA\Data\DynamicsResults\2019-12-10-4xVasamEGFPVK22Homo-HbNbHomo-Anterior-27_5C-NC14Peak';
load([DropboxFolder,filesep,'FrameInfo.mat'], 'FrameInfo')
nucleusDiameter = 3.96;
% Added by NL and GM on 11/23/2019
frameNumber =1;
xDim = FrameInfo(1).PixelsPerLine * FrameInfo(1).PixelSize;
yDim = FrameInfo(1).LinesPerFrame * FrameInfo(1).PixelSize;
if yDim > 150 && xDim > 150
    I = imread(names{frameNumber});
    pixelvalues = unique(I(:));
    threshold = pixelvalues(2);
    f_sigma = round(nucleusDiameter / FrameInfo(1).PixelSize);
    bwfill = imfill(I>threshold,'holes');
    I_inside = bwselect(bwfill,round(size(bwfill,2)/2), round(size(bwfill, 1)/2));
    I_inside=uint16((2^16-1)*I_inside);
    I_blurred = imfilter(I_inside,...
         fspecial('gaussian',2*blurSigma,blurSigma),'symmetric','conv');
    level = graythresh(I_blurred);
    I_mask = im2bw(I_blurred,level)
    embryoMask = imbinarize(I_blurred);    
else    
    if ~exist('embryoMask','var') || isempty(embryoMask)
        embryoMask = true(size(imread(names{frameNumber})));
    end
end

findNuclei(FrameInfo, names, frameNumber, nucleusDiameter, embryoMask)

%%  %% Get an embryo mask
% pixelvalues = unique(im2(:));
% 
% threshold = pixelvalues(2);
% blurSigma = round(6/PixelSize,0);
% I = im2;
% bwfill=imfill(I>threshold,'holes');
% I_inside = bwselect(bwfill,round(size(bwfill,2)/2),round(size(bwfill,1)/2));
% I_inside=uint16((2^16-1)*I_inside);
% I_blurred = imfilter(I_inside,...
%         fspecial('gaussian',2*blurSigma,blurSigma),'symmetric','conv');
% level = graythresh(I_blurred);
% I_mask = im2bw(I_blurred,level);
% figure(1)
% imagesc(I_mask)
% % im2 = imflatfield(I, blurSigma*4, I_mask);
% % figure(2)
% % imagesc(im2)
% 
% im2f = im2-imfilter(I,fspecial('disk', blurSigma/10));
% 
% 
% im2FF = imfilter(I,fspecial('disk',blurSigma*3));
% figure(3)
% imagesc(im2FF)
% 
% im2corrected = im2.*I_mask./im2FF;
% figure(4)
% imagesc(im2corrected)
% 
% 
% %% 
% load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo')
