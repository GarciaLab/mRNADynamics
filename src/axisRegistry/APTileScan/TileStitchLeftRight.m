%Hardcoded variables for initial implementation 
clear all, close all
SourcePath = 'E:/Gabriella/LivemRNA\Data\RawDynamicsData';
Prefix = '2019-12-09-4xVasamEGFPVK22Homo-HbNbHomo-Anterior-27_5C-NC13NC14';
Date = '2019-12-09';
EmbryoName = '4xVasamEGFPVK22Homo-HbNbHomo-Anterior-27_5C-NC13NC14';
SurfTileDataFolder = [SourcePath, filesep, Date, filesep, EmbryoName,filesep,...
    'FullEmbryo\SurfFileStacks'];

LeftFile = [SurfTileDataFolder,filesep,'SurfTile1.tif'];
RightFile = [SurfTileDataFolder,filesep,'SurfTile2.tif'];
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
TileIndex = 2;
T2ImageSlices = generateHisSlicesTileScan(LIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile2 = max(T2ImageSlices,[],3);
TileIndex = 3;
T3ImageSlices = generateHisSlicesTileScan(LIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile3 = max(T3ImageSlices, [], 3);
TileIndex = 4;
T4ImageSlices = generateHisSlicesTileScan(LIFImages, NSlices, NChannels, fiducialChannel, framesIndex, TileIndex);
tile4 = max(T4ImageSlices, [], 3);
xo1 = [];
yo1 = [];
% x and y limits
x_min=0; x_max=100;
y_min=-100; y_max=100;



% x and y limits
x_min=0; x_max=150;
y_min=-150; y_max=150;




    
[xo12 yo12] = autoStitchLR(tile1, tile2, x_max, y_max, x_min, y_min);
imshow(imstitchLR(tile1, tile2, xo12, yo12, [1 2]));




[xo34 yo34] = autoStitchLR(tile3, tile4, x_max, y_max, x_min, y_min);
figure(2)
imshow(imstitchLR(tile3, tile4, xo34, yo34, [1 2]));


% x and y limits
x_min=-150; x_max=150;
y_min=0; y_max=150;

[xo13 yo13] = autoStitchBT(tile1, tile3, x_max, y_max, x_min, y_min);
figure(3)
imshow(imstitchBT(tile1, tile3, xo13, yo13, [1 2]));


[xo24 yo24] = autoStitchBT(tile2, tile4, x_max, y_max, x_min, y_min);
figure(4)
imshow(imstitchBT(tile2, tile4, xo24, yo24, [1 2]));





%% 



function [xo1 yo1] = autoStitchLR(left, right, x_max, y_max, x_min, y_min)
    if x_min < 0
        error(['x offset cannot be negative: tiled images must have an overlap',...
        'greater than or equal to zero. Choose a non-negative value for x_min.' ]);
    end
    XRange = x_min:1:x_max;
    YRange = y_min:1:y_max;

    leftN = normalizer(left);
    rightN = normalizer(right);

    [xo1, yo1] = matchManualLR(leftN, rightN, YRange, XRange);
end

function [xo1 yo1] = autoStitchBT(bottom, top, x_max, y_max, x_min, y_min)
    if y_min < 0
        error(['y offset cannot be negative: tiled images must have an overlap',...
        'greater than or equal to zero. Choose a non-negative value for y_min.' ]);
    end
    XRange = x_min:1:x_max;
    YRange = y_min:1:y_max;

    topN = normalizer(top);
    bottomN = normalizer(bottom);

    [xo1, yo1] = matchManualBT(bottomN, topN, YRange, XRange);
end


%Works
function [xoffset, yoffset] = matchManualLR(leftN, rightN, yRange, xRange) 
    % goes through where the widths are close, and match the two normalized
    % images to see how good they match
    ys = repmat(yRange, 1, length(xRange));
    xs = reshape(repmat(xRange, length(yRange), 1), 1, []);
    
    % claculate scores in region and find the best
    scores = arrayfun(@(x, y) picDiffLR(leftN, rightN, x, y), xs, ys);
    %areas = arrayfun(@(x, y) calcAreaLR(x, y, size(leftN, 2), size(leftN, 1)), xs, ys);
    best = find(scores==min(scores), 1, 'first');
    
    xoffset = xs(best);
    yoffset = ys(best);
end



function [xoffset, yoffset] = matchManualBT(bottomN, topN, yRange, xRange) 
    % goes through where the widths are close, and match the two normalized
    % images to see how good they match

    ys = repmat(yRange, 1, length(xRange));
    xs = reshape(repmat(xRange, length(yRange), 1), 1, []);
    
    % claculate scores in region and find the best
    scores = arrayfun(@(x, y) picDiffBT(bottomN, topN, x, y), xs, ys);
    best = find(scores==min(scores), 1, 'first');
    
    xoffset = xs(best);
    yoffset = ys(best);
end

function overlap = calcAreaLR(x, y, w, h)
    overlap = (h-abs(y))*x;
end

function overlap = calcAreaBT(x, y, w, h)
    overlap = y*(w-abs(x));
end

%Works
function score = picDiffLR(leftN, rightN, x, y)
    % returns how well the two normalized images match.
    [h, w] = size(leftN);
    if y>=0
        leftWindow = leftN(y+1:h,w-x+1:w);
        rightWindow = rightN(1:h-y, 1:x);
        area = (h-y)*x;
    else
        leftWindow = leftN(1:h+y,w-x+1:w);
        rightWindow = rightN(-y+1:h,1:x);
        area = (h+y)*x;
    end
    
%         if x>=0,
%         leftWindow = leftN(h-y+1:h, x+1:w);
%         rightWindow = rightN(1:y, 1:w-x);
%         area = (w-x)*y;
%     else
%         leftWindow = leftN(h-y+1:h, 1:w+x);
%         rightWindow = rightN(1:y, -x+1:w);
%         area = (w+x)*y;
%     end
    
    diffSquared = (leftWindow - rightWindow) .^ 2;
    score = sqrt(sum(diffSquared(:)))/area;
end

%Works
function score = picDiffBT(bottomN, topN, x, y)
    % returns how well the two normalized images match.
    [h, w] = size(topN);
    % 
    if x>=0 % Case where bottom is shifted to the right relative to top 
        topWindow = topN(h-y+1:h,x+1:w);
        bottomWindow = bottomN(1:y,1:w-x);
        area = (w-x)*y;
        %leftWindow = leftN(y+1:h,w-x+1:w);
        %rightWindow = rightN(1:h-y, 1:x);
        %area = (h-y)*x;
    else
        topWindow = topN(h-y+1:h,1:w+x);
        bottomWindow = bottomN(1:y, -x+1:w);
        area = (w+x)*y;
        %leftWindow = leftN(1:h+y,w-x+1:w);
        %rightWindow = rightN(-y+1:h,1:x);
        %area = (h+y)*x;
    end
    
    diffSquared = (bottomWindow - topWindow) .^ 2;
    score = sum(diffSquared(:))/area;
end

%Works
function im2 = normalizer(im, width) 
    if ~exist('width', 'var'), width = 80; end
    im = double(im);
    im = im - imfilter(im, fspecial('average', width)); %mean normalize
    im = im ./ (sqrt(imfilter(im .* im, fspecial('average', width))));
    im = im / std(im(:));
    im2 = double(im); % > mean(im(:)) + 2*std(im(:)));
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
