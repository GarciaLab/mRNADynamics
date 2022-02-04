function [FrameIndex, SpotIndex,NumFrames, NumGoodSpots] = Store3DNanocageInfo(Prefix, MaxPixelValue, FrameIndex, SpotIndex, UseFilteredImages)
if ~exist('MaxPixelValue', 'var')
    MaxPixelValue = 20;
end
if ~exist('UseFilteredImages', 'var')
    UseFilteredImages = true;
end
liveExperiment = LiveExperiment(Prefix);
Spots = getSpots(liveExperiment);
snippet_size = ceil(1300/(2*liveExperiment.pixelSize_nm));
max_snippet_size = snippet_size*2;
zoomed_out_snippet_size = snippet_size*10;
xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
NumFrames = length(Spots);
TifList = dir([liveExperiment.preFolder, filesep, '*ch01.tif']);
ProbList = dir([liveExperiment.procFolder,'dogs',  filesep, 'prob*']);
MovieCells = cell(1, NumFrames);
ProbCells = cell(1, NumFrames);
SliceCells = cell(1,NumFrames);
NumSlices = 0;
for i = 1:NumFrames
    MovieCells{i} = imreadStack([TifList(i).folder, filesep,  TifList(i).name]);
    SliceCells{i} = imreadStack([TifList(i).folder, filesep,  TifList(i).name]);
    ProbCells{i} = imreadStack([ProbList(i).folder, filesep,  ProbList(i).name]);
    NumSlices = NumSlices + size(MovieCells{i}, 3);
end


%%

% AllSpotPoints = cell(1,NumFrames);
% for FrameIndex = 1:NumFrames
%     AllSpotPoints{FrameIndex}.rows = [];
%     AllSpotPoints{FrameIndex}.cols = [];
%     AllSpotPoints{FrameIndex}.z = [];
%     for SpotIndex = 1:length(Spots(FrameIndex).Fits)
%         for zIndex = 1:length(Spots(FrameIndex).Fits(SpotIndex).z)
%             CurrentZ = Spots(FrameIndex).Fits(SpotIndex).z(zIndex);
%             PixelList = Spots(FrameIndex).Fits(SpotIndex).PixelList{zIndex};
%             AllSpotPoints{FrameIndex}.rows =[AllSpotPoints{FrameIndex}.rows PixelList(:,1).'];
%             AllSpotPoints{FrameIndex}.cols =[AllSpotPoints{FrameIndex}.cols PixelList(:,2).'];
%             AllSpotPoints{FrameIndex}.z =[AllSpotPoints{FrameIndex}.z CurrentZ*ones(1, size(PixelList, 1), 'uint8')];
%         end
%
%
%     end
% end



%%

ValidFrameIndices = [];
for fr = 1:length(Spots)
    if ~isempty(Spots(fr).Fits)
        ValidFrameIndices = [ValidFrameIndices fr];
    end
end



if isfile([liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat'])
    load([liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat']);
    if ~exist('BestIndices', 'var')
        BestIndices = cell(1, length(Spots));
    end
    if ~exist('AggregateIndices', 'var')
        AggregateIndices = cell(1, length(Spots));
    end
    if ~exist('Cat1Indices', 'var')
        Cat1Indices = cell(1, length(Spots));
    end
    if ~exist('Cat2Indices', 'var')
        Cat2Indices = cell(1, length(Spots));
     end
    if ~exist('Cat3Indices', 'var')
        Cat3Indices = cell(1, length(Spots));
    end
    if ~exist('Cat4Indices', 'var')
        Cat4Indices = cell(1, length(Spots));
    end
    if ~exist('Cat5Indices', 'var')
        Cat5Indices = cell(1, length(Spots));
    end
else
    GoodIndices = cell(1, length(Spots));
    BadIndices = cell(1, length(Spots));
    MultiIndices = cell(1, length(Spots));
    BestIndices = cell(1, length(Spots));
    AggregateIndices = cell(1, length(Spots));
    Cat1Indices = cell(1, length(Spots));
    Cat2Indices = cell(1, length(Spots));
    Cat3Indices = cell(1, length(Spots));
    Cat4Indices = cell(1, length(Spots));
    Cat5Indices = cell(1, length(Spots));
end

MaxStackSize = 0;
for fr = 1:length(Spots)
    MaxStackSize = max([MaxStackSize, max([Spots(fr).Fits(:).zCount])]);
end
MaxStackSize = min([MaxStackSize, 8]);


SubplotPositionInfo = GetSubplotPositioningParameters(MaxStackSize, false);
if ~exist('FrameIndex', 'var')
    FrameIndex = ValidFrameIndices(1);
end
if ~exist('SpotIndex', 'var')
    SpotIndex = 1;
end

%%
mean_rho = 0;
mean_sigx = 0.5;
mean_sigy = 0.5;
filter_size = 5;
middle = 3;
custom_filter = zeros(5,5,'double');
for i = 1:filter_size
    for j = 1:filter_size
        custom_filter(i,j) = exp(-1/(2*(1-mean_rho)^2)*(((j-middle)^2)/(2*mean_sigx^2)+...
            ((i-middle)^2)/(2*mean_sigy^2)-(2*mean_rho*(j-middle)*(i-middle))/(mean_sigx*mean_sigy)));
    end
end


if UseFilteredImages
    for i = 1:NumFrames
        for j =2:(size(MovieCells{i},3)-1)
            SliceCells{i}(:,:,j) = imfilter(SliceCells{i}(:,:,j),custom_filter);
        end
    end
end
%%
theta = 0 : (2 * pi / 10000) : (2 * pi);
pline_x = snippet_size/2 * cos(theta) + snippet_size+1;
pline_y = snippet_size/2  * sin(theta) + snippet_size+1;
max_pline_x = snippet_size/2 * cos(theta) + max_snippet_size+1;
max_pline_y = snippet_size/2  * sin(theta) + max_snippet_size+1;
zoomed_pline_x = snippet_size/2 * cos(theta) + zoomed_out_snippet_size+1;
zoomed_pline_y = snippet_size/2  * sin(theta) + zoomed_out_snippet_size+1;
close all
FigAx = cell(1, MaxStackSize);

TempFigure = figure(1);

set(TempFigure,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
set(gcf,'color','w');
SnippetImage2 = zeros(2*snippet_size+1,2*snippet_size+1, MaxStackSize, 'double');
SnippetImage3 = zeros(2*max_snippet_size+1,2*max_snippet_size+1, MaxStackSize, 'double');
SnippetImage4 = zeros(2*zoomed_out_snippet_size+1,2*zoomed_out_snippet_size+1, MaxStackSize, 'double');
clim = [0, MaxPixelValue];

for SubplotIndex = 1:MaxStackSize
    FigAx{SubplotIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotIndex);
    
    imagesc(FigAx{SubplotIndex} , squeeze(SnippetImage2(:,:,1)), clim);
    if SubplotIndex < MaxStackSize
        colormap(FigAx{SubplotIndex},'gray')
    else
        colormap(FigAx{SubplotIndex},'gray')%,'parula')
    end
    axis off
end

for SubplotIndex =1:MaxStackSize
    pos = get(FigAx{SubplotIndex} , 'position');
    pos(1) = SubplotPositionInfo.SubplotXPositions(SubplotIndex);
    pos(3) = SubplotPositionInfo.SubplotWidth;
    pos(2) = SubplotPositionInfo.SubplotYPositions(SubplotIndex);
    pos(4) = SubplotPositionInfo.SubplotHeight;
    set(FigAx{SubplotIndex}, 'position', pos);
end
CurrentZ = Spots(FrameIndex).Fits(SpotIndex).z(1);
plot_title= {strrep(TifList(FrameIndex).name, '_', '\_'), ['Spot Index : ', num2str(SpotIndex)]};

MiddleZ = uint8(floor(median(Spots(FrameIndex).Fits(SpotIndex).z)));
RefIndex = find(Spots(FrameIndex).Fits(SpotIndex).z == MiddleZ);
%%
xmin = Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex)-snippet_size;
xmax = Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex)+snippet_size;
ymin = Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex)-snippet_size;
ymax = Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex)+snippet_size;
if xmin> 0 & ymin > 0 & xmax <= xDim & ymax <= yDim
    SnippetImage2(:,:, 1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = ...
        SliceCells{FrameIndex}(ymin:ymax,xmin:xmax, Spots(FrameIndex).Fits(SpotIndex).z);
else
    if xmin <= 0
        xSnip = 1-(xmin);
    else
        xSnip =1;
    end
    if ymin <= 0
        ySnip = 1-(ymin);
    else
        ySnip=1;
    end
    xmaxSnip = xSnip + (xmax-xmin);
    ymaxSnip = ySnip + (ymax-ymin);
    SnippetImage2(ySnip:ymaxSnip, xSnip:xmaxSnip, 1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = SliceCells{FrameIndex}(ymin:ymax,xmin:xmax, Spots(FrameIndex).Fits(SpotIndex).z);
end

zoomed_xmin = double(Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex))-zoomed_out_snippet_size;
zoomed_xmax = double(Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex))+zoomed_out_snippet_size;
zoomed_ymin = double(Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex))-zoomed_out_snippet_size;
zoomed_ymax = double(Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex))+zoomed_out_snippet_size;
if zoomed_xmin> 0 & zoomed_ymin > 0 & zoomed_xmax <= xDim & zoomed_ymax <= yDim
    SnippetImage4(:,:, 1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = ...
        MovieCells{FrameIndex}(zoomed_ymin:zoomed_ymax,zoomed_xmin:zoomed_xmax, Spots(FrameIndex).Fits(SpotIndex).z);
else
    if zoomed_xmin <= 0
        zoomed_xSnip = 1-(zoomed_xmin);
        zoomed_xmin = 1;
    else
        zoomed_xSnip =1;
    end
    if zoomed_ymin <= 0
        zoomed_ySnip = 1-(zoomed_ymin);
        zoomed_ymin = 1;
    else
        zoomed_ySnip=1;
    end
    
    
    if zoomed_xmax > xDim
        zoomed_xmax = xDim;
    end
    if zoomed_ymax > yDim
        zoomed_ymax = yDim;
    end
    zoomed_xmaxSnip = zoomed_xSnip + (zoomed_xmax-zoomed_xmin);
    zoomed_ymaxSnip = zoomed_ySnip + (zoomed_ymax-zoomed_ymin);
    SnippetImage4(zoomed_ySnip:zoomed_ymaxSnip, zoomed_xSnip:zoomed_xmaxSnip, 1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = MovieCells{FrameIndex}(zoomed_ymin:zoomed_ymax,zoomed_xmin:zoomed_xmax, Spots(FrameIndex).Fits(SpotIndex).z);
end


max_xmin = double(Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex))-max_snippet_size;
max_xmax = double(Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex))+max_snippet_size;
max_ymin = double(Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex))-max_snippet_size;
max_ymax = double(Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex))+max_snippet_size;
if max_xmin> 0 & max_ymin > 0 & max_xmax <= xDim & max_ymax <= yDim
    SnippetImage3(:,:, 1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = ...
        MovieCells{FrameIndex}(max_ymin:max_ymax,max_xmin:max_xmax, Spots(FrameIndex).Fits(SpotIndex).z);
else
    if max_xmin <= 0
        max_xSnip = 1-(max_xmin);
    else
        max_xSnip =1;
    end
    if max_ymin <= 0
        max_ySnip = 1-(max_ymin);
    else
        max_ySnip=1;
    end
    
    if max_xmax > xDim
        max_xmax = xDim;
    end
    if max_ymax > yDim
        max_ymax = yDim;
    end
    
    
    max_xmaxSnip = max_xSnip + (max_xmax-max_xmin);
    max_ymaxSnip = max_ySnip + (max_ymax-max_ymin);
    SnippetImage3(max_ySnip:max_ymaxSnip, max_xSnip:max_xmaxSnip, 1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = MovieCells{FrameIndex}(max_ymin:max_ymax,max_xmin:max_xmax, Spots(FrameIndex).Fits(SpotIndex).z);
end


ScatterHandlesRed =cell(1,MaxStackSize);
ScatterHandlesYellow =cell(1,MaxStackSize);
ScatterHandlesBlue =cell(1,MaxStackSize);
for SubplotIndex = 1:MaxStackSize-2
    set(TempFigure, 'currentaxes', FigAx{SubplotIndex});
    set( FigAx{SubplotIndex}, 'YLimMode', 'manual' )
    set( FigAx{SubplotIndex}, 'XLimMode', 'manual' )
    xlim([0.5 2*snippet_size+1.5]);
    ylim([0.5 2*snippet_size+1.5]);
    if SubplotIndex <= length(Spots(FrameIndex).Fits(SpotIndex).z)
        
        imagesc(FigAx{SubplotIndex} , squeeze(SnippetImage2(:,:,SubplotIndex)), clim);
        title(['z = ', num2str(Spots(FrameIndex).Fits(SpotIndex).z(SubplotIndex))])
        hold on
        colormap(FigAx{SubplotIndex},'gray')
        axis off
        
        plot(FigAx{SubplotIndex}, pline_x, pline_y, '.r');
        row_pixels= Spots(FrameIndex).Fits(SpotIndex).PixelList{SubplotIndex}(:,1).'-double(xmin)+1;
        col_pixels= Spots(FrameIndex).Fits(SpotIndex).PixelList{SubplotIndex}(:,2).'-double(ymin)+1;
        
        
        
        
        ScatterHandlesRed{SubplotIndex} = scatter(row_pixels,col_pixels, 'r*', 'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
        [mid_row, mid_cols] = find(squeeze(SnippetImage2(:,:,SubplotIndex)) >= MaxPixelValue/2);
        ScatterHandlesYellow{SubplotIndex} = scatter(mid_cols,mid_row,  'y*',  'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
        [high_row, high_cols] = find(squeeze(SnippetImage2(:,:,SubplotIndex)) >= MaxPixelValue);
        ScatterHandlesBlue{SubplotIndex} = scatter(high_cols,high_row,  'c*',  'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
        
        hold off
    else
        imagesc(FigAx{SubplotIndex}, zeros(2*snippet_size+1,2*snippet_size+1,'double'), clim);
        hold on
        plot(FigAx{SubplotIndex}, pline_x, pline_y, '.r');
        [high_row, high_cols] = find(squeeze(SnippetImage2(:,:,SubplotIndex)) > MaxPixelValue);
        [mid_row, mid_cols] = find(squeeze(SnippetImage2(:,:,SubplotIndex)) > MaxPixelValue/2);
        ScatterHandlesRed{SubplotIndex} = scatter(high_cols,high_row, 'r*', 'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
        
        ScatterHandlesYellow{SubplotIndex} = scatter(mid_cols, mid_row,  'y*',  'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
        ScatterHandlesBlue{SubplotIndex} = scatter(high_cols,high_row,  'c*',  'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
        
        colormap(FigAx{SubplotIndex},'gray')
        axis off
        hold off
        title('')
    end
    
    
    
end

set(TempFigure, 'currentaxes', FigAx{MaxStackSize-1});
set( FigAx{MaxStackSize-1}, 'YLimMode', 'manual' )
set( FigAx{MaxStackSize-1}, 'XLimMode', 'manual' )
xlim([0.5 2*max_snippet_size+1.5]);
ylim([0.5 2*max_snippet_size+1.5]);
imagesc(FigAx{MaxStackSize-1} , sum(SnippetImage3,3), clim*1);% max(SnippetImage3,[],3), clim);
hold on

title('Max Proj.')
colormap(FigAx{MaxStackSize-1},'gray')%,'parula')
axis off
plot(FigAx{MaxStackSize-1}, max_pline_x, max_pline_y, '.r');
hold off

set(TempFigure, 'currentaxes', FigAx{MaxStackSize});
set( FigAx{MaxStackSize}, 'YLimMode', 'manual' )
set( FigAx{MaxStackSize}, 'XLimMode', 'manual' )
xlim([0.5 2*zoomed_out_snippet_size+1.5]);
ylim([0.5 2*zoomed_out_snippet_size+1.5]);
imagesc(FigAx{MaxStackSize} , sum(SnippetImage4,3), clim*1);% max(SnippetImage3,[],3), clim);
hold on

title('Max Proj.')
colormap(FigAx{MaxStackSize},'parula')%,'parula')
axis off
plot(FigAx{MaxStackSize}, zoomed_pline_x, zoomed_pline_y,'.r',...
    'MarkerSize',3);
rowVector = Spots(FrameIndex).Fits(SpotIndex).xDoG-double(zoomed_xmin)+1;
colVector = Spots(FrameIndex).Fits(SpotIndex).yDoG-double(zoomed_ymin)+1;


%ScatterHandlesRed{MaxStackSize} = scatter(FigAx{MaxStackSize}, rowVector, colVector,10, 'MarkerFaceColor','r','MarkerEdgeColor','r');
hold off


%%
sgtitle(plot_title)

if isempty(BadIndices{FrameIndex})
    BadIndices{FrameIndex} = [];
end

if isempty(GoodIndices{FrameIndex})
    GoodIndices{FrameIndex} = [];
end

if isempty(MultiIndices{FrameIndex})
    MultiIndices{FrameIndex} = [];
end

if isempty(BestIndices{FrameIndex})
    BestIndices{FrameIndex} = [];
end

if isempty(AggregateIndices{FrameIndex})
    AggregateIndices{FrameIndex} = [];
end
if isempty(Cat1Indices{FrameIndex})
            Cat1Indices{FrameIndex} = [];
        end
        if isempty(Cat2Indices{FrameIndex})
            Cat2Indices{FrameIndex} = [];
        end
        if isempty(Cat3Indices{FrameIndex})
            Cat3Indices{FrameIndex} = [];
        end
        if isempty(Cat4Indices{FrameIndex})
            Cat4Indices{FrameIndex} = [];
        end
        if isempty(Cat5Indices{FrameIndex})
            Cat5Indices{FrameIndex} = [];
        end
        

if ismember(SpotIndex, BestIndices{FrameIndex})
    set(TempFigure,'Color','c')
elseif ismember(SpotIndex, AggregateIndices{FrameIndex})
    set(TempFigure,'Color','m')
elseif ismember(SpotIndex, BadIndices{FrameIndex})
    set(TempFigure,'Color','r')
elseif ismember(SpotIndex, GoodIndices{FrameIndex})
    set(TempFigure,'Color','g')
elseif ismember(SpotIndex, MultiIndices{FrameIndex})
    set(TempFigure,'Color','b')
else
    set(TempFigure,'Color','y')
end

hold off

%%

currentCharacter = 0;
set(0, 'CurrentFigure', TempFigure)


while (currentCharacter~='x')
    CurrentZ = Spots(FrameIndex).Fits(SpotIndex).z(1);
    NumBestSpots = length(BestIndices{FrameIndex});
    NumGoodSpots = length(GoodIndices{FrameIndex});
    NumBadSpots = length(BadIndices{FrameIndex});
    NumMultiSpots = length(MultiIndices{FrameIndex});
    NumAggregates = length(AggregateIndices{FrameIndex});
    NumCat1s = length(Cat1Indices{FrameIndex});
    NumCat2s = length(Cat2Indices{FrameIndex});
    NumCat3s = length(Cat3Indices{FrameIndex});
    NumCat4s = length(Cat4Indices{FrameIndex});
    NumCat5s = length(Cat5Indices{FrameIndex});
    SpotCat = 0;
    if ismember(SpotIndex, Cat1Indices{FrameIndex})
        SpotCat = 1;
    elseif ismember(SpotIndex, Cat2Indices{FrameIndex})
        SpotCat = 2;
    elseif ismember(SpotIndex, Cat3Indices{FrameIndex})
        SpotCat = 3;
    elseif ismember(SpotIndex, Cat4Indices{FrameIndex})
        SpotCat = 4;
    elseif ismember(SpotIndex, Cat5Indices{FrameIndex})
        SpotCat = 5;
    end
    plot_title= {[strrep(TifList(FrameIndex).name, '_', '\_'), ' (Frame Count: ', num2str(length(Spots)), ')' ], ['Spot Index: ',...
        num2str(SpotIndex),'/',num2str(length(Spots(FrameIndex).Fits)),...
        ', Spot Cat.: ', num2str(SpotCat),...
        ', Z Count: ', num2str(length(Spots(FrameIndex).Fits(SpotIndex).z)),...
        ', Good: ', num2str(NumGoodSpots), ', Bad: ',...
        num2str(NumBadSpots), ', Multi: ', num2str(NumMultiSpots),...
        ', Best: ', num2str(NumBestSpots),...
        ', Agg: ', num2str(NumAggregates),...
        ', Cat1: ', num2str(NumCat1s),...
        ', Cat2: ', num2str(NumCat2s),...
        ', Cat3: ', num2str(NumCat3s),...
        ', Cat4: ', num2str(NumCat4s),...
        ', Cat5: ', num2str(NumCat5s) ...
        ]};
    
    MiddleZ = uint8(floor(median(Spots(FrameIndex).Fits(SpotIndex).z)));
    RefIndex = find(Spots(FrameIndex).Fits(SpotIndex).z == MiddleZ);
    %Load subsequent images
    
    
    
    
    SnippetImage = zeros(2*snippet_size+1,2*snippet_size+1,MaxStackSize, 'double');
    xmin = Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex)-snippet_size;
    xmax = Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex)+snippet_size;
    ymin = Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex)-snippet_size;
    ymax = Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex)+snippet_size;
    if xmin> 0 & ymin > 0 & xmax <= xDim & ymax <= yDim
        SnippetImage = SliceCells{FrameIndex}(ymin:ymax,xmin:xmax,Spots(FrameIndex).Fits(SpotIndex).z);%.*double(SpotMat(xmin:xmax,ymin:ymax,ReorderedIndex));
    else
        if xmin <= 0
            xSnip = 1-(xmin);
            xmin = 1;
        else
            xSnip =1;
        end
        if ymin <= 0
            ySnip = 1-(ymin);
            ymin = 1;
        else
            ySnip=1;
        end
        if xmax > xDim
            xmax = xDim;
        end
        if ymax > yDim
            ymax = yDim;
        end
        xmaxSnip = xSnip + (xmax-xmin);
        ymaxSnip = ySnip + (ymax-ymin);
        SnippetImage(ySnip:ymaxSnip,xSnip:xmaxSnip,  1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = SliceCells{FrameIndex}(ymin:ymax,xmin:xmax,Spots(FrameIndex).Fits(SpotIndex).z);%.*double(SpotMat(xmin:xmax,ymin:ymax,ReorderedIndex));
    end
    
    SnippetImage3 = zeros(2*max_snippet_size+1,2*max_snippet_size+1,MaxStackSize, 'double');
    
    max_xmin = double(Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex))-max_snippet_size;
    max_xmax = double(Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex))+max_snippet_size;
    max_ymin = double(Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex))-max_snippet_size;
    max_ymax = double(Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex))+max_snippet_size;
    if max_xmin> 0 & max_ymin > 0 & max_xmax <= xDim & max_ymax <= yDim
        SnippetImage3 = MovieCells{FrameIndex}(max_ymin:max_ymax,max_xmin:max_xmax,Spots(FrameIndex).Fits(SpotIndex).z);%.*double(SpotMat(xmin:xmax,ymin:ymax,ReorderedIndex));
    else
        if max_xmin <= 0
            max_xSnip = 1-(max_xmin);
            max_xmin = 1;
        else
            max_xSnip =1;
        end
        if max_ymin <= 0
            max_ySnip = 1-(max_ymin);
            max_ymin = 1;
        else
            max_ySnip=1;
        end
        if max_xmax > xDim
            max_xmax = xDim;
        end
        if max_ymax > yDim
            max_ymax = yDim;
        end
        max_xmaxSnip = max_xSnip + (max_xmax-max_xmin);
        max_ymaxSnip = max_ySnip + (max_ymax-max_ymin);
        SnippetImage3(max_ySnip:max_ymaxSnip,max_xSnip:max_xmaxSnip,  1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = MovieCells{FrameIndex}(max_ymin:max_ymax,max_xmin:max_xmax,Spots(FrameIndex).Fits(SpotIndex).z);%.*double(SpotMat(xmin:xmax,ymin:ymax,ReorderedIndex));
    end
    
    SnippetImage4 = zeros(2*zoomed_out_snippet_size+1,2*zoomed_out_snippet_size+1,MaxStackSize, 'double');
    zoomed_xmin = double(Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex))-zoomed_out_snippet_size;
    zoomed_xmax = double(Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex))+zoomed_out_snippet_size;
    zoomed_ymin = double(Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex))-zoomed_out_snippet_size;
    zoomed_ymax = double(Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex))+zoomed_out_snippet_size;
    if zoomed_xmin> 0 & zoomed_ymin > 0 & zoomed_xmax <= xDim & zoomed_ymax <= yDim
        SnippetImage4(:,:, 1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = ...
            MovieCells{FrameIndex}(uint16(zoomed_ymin):uint16(zoomed_ymax),uint16(zoomed_xmin):uint16(zoomed_xmax), Spots(FrameIndex).Fits(SpotIndex).z);
    else
        if zoomed_xmin <= 0
            zoomed_xSnip = 1-(zoomed_xmin);
            zoomed_xmin = 1;
        else
            zoomed_xSnip =1;
        end
        if zoomed_ymin <= 0
            zoomed_ySnip = 1-(zoomed_ymin);
            zoomed_ymin = 1;
        else
            zoomed_ySnip=1;
        end
        
        if zoomed_xmax > xDim
            zoomed_xmax = xDim;
        end
        if zoomed_ymax > yDim
            zoomed_ymax = yDim;
        end
        
        zoomed_xmaxSnip = zoomed_xSnip + (zoomed_xmax-zoomed_xmin);
        zoomed_ymaxSnip = zoomed_ySnip + (zoomed_ymax-zoomed_ymin);
        SnippetImage4(zoomed_ySnip:zoomed_ymaxSnip, zoomed_xSnip:zoomed_xmaxSnip, 1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = MovieCells{FrameIndex}(zoomed_ymin:zoomed_ymax,zoomed_xmin:zoomed_xmax, Spots(FrameIndex).Fits(SpotIndex).z);
    end
    
    
    for SubplotIndex = 1:MaxStackSize-2
        set(TempFigure, 'currentaxes', FigAx{SubplotIndex});
        %hold(FigAx{SubplotIndex}, 'on')
        if SubplotIndex <= length(Spots(FrameIndex).Fits(SpotIndex).z)
            FigAx{SubplotIndex}.Children(end).CData = squeeze(SnippetImage(:,:,SubplotIndex));
            FigAx{SubplotIndex}.Title.String = ['z = ', num2str(Spots(FrameIndex).Fits(SpotIndex).z(SubplotIndex))];
            %             CurrentZ = uint8(Spots(FrameIndex).Fits(SpotIndex).z(SubplotIndex));
            %             allowed_spots = (AllSpotPoints{FrameIndex}.z == CurrentZ) & ...
            %                 (AllSpotPoints{FrameIndex}.rows >= xmin) & ...
            %                 (AllSpotPoints{FrameIndex}.rows <= xmax) &...
            %                 (AllSpotPoints{FrameIndex}.cols >= ymin) & ...
            %                 (AllSpotPoints{FrameIndex}.cols <= ymax);
            %             row_pixels2 = AllSpotPoints{FrameIndex}.rows(allowed_spots)-double(xmin)+1;
            %             col_pixels2 = AllSpotPoints{FrameIndex}.cols(allowed_spots)-double(ymin)+1;
            %             scatter(row_pixels2,col_pixels2, 'c*', 'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
            row_pixels= Spots(FrameIndex).Fits(SpotIndex).PixelList{SubplotIndex}(:,1).'-double(xmin)+1;
            col_pixels= Spots(FrameIndex).Fits(SpotIndex).PixelList{SubplotIndex}(:,2).'-double(ymin)+1;
            set(ScatterHandlesRed{SubplotIndex}, 'XData', row_pixels, 'YData', col_pixels);
            set(FigAx{SubplotIndex}.Children(end-2),'Visible','on'); %'off' or 'on'
            [mid_row, mid_cols] = find(squeeze(SnippetImage(:,:,SubplotIndex)) >= MaxPixelValue/2);
            [high_row, high_cols] = find(squeeze(SnippetImage(:,:,SubplotIndex)) >= MaxPixelValue);
            set(ScatterHandlesYellow{SubplotIndex},'XData', mid_cols,'YData', mid_row);
            set(ScatterHandlesBlue{SubplotIndex},'XData', high_cols,'YData', high_row);
            set(FigAx{SubplotIndex}.Children(end-3),'Visible','on'); %'off' or 'on'
            set(FigAx{SubplotIndex}.Children(end-4),'Visible','on'); %'off' or 'on'
            hold off
            
        else
            FigAx{SubplotIndex}.Children(end).CData = zeros(2*snippet_size+1,2*snippet_size+1,'double');
            FigAx{SubplotIndex}.Title.String = '';
            
            %[high_row, high_cols] = find(FigAx{SubplotIndex}.Children(end).CData  > MaxPixelValue);
            %set(ScatterHandlesRed{SubplotIndex},'XData', high_cols,'YData', high_row);
            %set(ScatterHandlesYellow{SubplotIndex},'XData', high_cols,'YData', high_row);
            set(FigAx{SubplotIndex}.Children(end-2),'Visible','off'); %'off' or 'on'
            set(FigAx{SubplotIndex}.Children(end-3),'Visible','off'); %'off' or 'on'
            set(FigAx{SubplotIndex}.Children(end-4),'Visible','off'); %'off' or 'on'
            hold off
        end
        %hold(FigAx{SubplotIndex}, 'off')
        
        
    end
    
    set(TempFigure, 'currentaxes', FigAx{MaxStackSize-1});
    FigAx{MaxStackSize-1}.Children(end).CData = sum(SnippetImage3,3);%max(SnippetImage3, [], 3);
    
    set(TempFigure, 'currentaxes', FigAx{MaxStackSize});
    FigAx{MaxStackSize}.Children(end).CData = sum(SnippetImage4,3);%max(SnippetImage3, [], 3);
    rowVector = Spots(FrameIndex).Fits(SpotIndex).xDoG-double(zoomed_xmin)+1;
    colVector = Spots(FrameIndex).Fits(SpotIndex).yDoG-double(zoomed_ymin)+1;
    
    
    % set(ScatterHandlesRed{MaxStackSize}, 'XData', rowVector, 'YData', colVector);
    %     plot(FigAx{MaxStackSize}, zoomed_pline_x, zoomed_pline_y,'.r',...
    %     'MarkerSize',5);
    hold off
    sgtitle(plot_title)
    if ismember(SpotIndex, BestIndices{FrameIndex})
        set(TempFigure,'Color','c')
    elseif ismember(SpotIndex, AggregateIndices{FrameIndex})
        set(TempFigure,'Color','m')
    elseif ismember(SpotIndex, BadIndices{FrameIndex})
        set(TempFigure,'Color','r')
    elseif ismember(SpotIndex, GoodIndices{FrameIndex})
        set(TempFigure,'Color','g')
    elseif ismember(SpotIndex, MultiIndices{FrameIndex})
        set(TempFigure,'Color','b')
        
    else
        set(TempFigure,'Color','y')
    end
    
    hold off
    %%
    
    ct=waitforbuttonpress;
    currentCharacter=get(TempFigure,'currentcharacter');
    currentMouse=get(TempFigure,'CurrentPoint');
    
    if (ct~=0)&(currentCharacter=='.')&(SpotIndex<length(Spots(FrameIndex).Fits))
        SpotIndex=SpotIndex+1;
    elseif (ct~=0)&(currentCharacter==',')&(SpotIndex>1)
        SpotIndex=SpotIndex-1;
    elseif (ct~=0)&(currentCharacter=='>')
        if ~isempty(GoodIndices{FrameIndex})
            if SpotIndex<max(GoodIndices{FrameIndex})
                SpotIndex = GoodIndices{FrameIndex}(find(GoodIndices{FrameIndex} > SpotIndex, 1));
            end
        end
     elseif (ct~=0)&(currentCharacter=='<')
        if ~isempty(GoodIndices{FrameIndex})
            if SpotIndex>min(GoodIndices{FrameIndex})
                SpotIndex = GoodIndices{FrameIndex}(find(GoodIndices{FrameIndex} < SpotIndex, 1, 'last'));
            end
        end     
    elseif (ct~=0)&(currentCharacter=='m')&(FrameIndex<max(ValidFrameIndices))
        FrameIndex=ValidFrameIndices(find(ValidFrameIndices > FrameIndex, 1));%FrameIndex+1;
        if isempty(BadIndices{FrameIndex})
            BadIndices{FrameIndex} = [];
        end
        
        if isempty(GoodIndices{FrameIndex})
            GoodIndices{FrameIndex} = [];
        end
        
        if isempty(MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = [];
        end
        
        if isempty(BestIndices{FrameIndex})
            BestIndices{FrameIndex} = [];
        end
        
        if isempty(AggregateIndices{FrameIndex})
            AggregateIndices{FrameIndex} = [];
        end
        SpotIndex = 1;
        
    elseif (ct~=0)&(currentCharacter=='z')&(FrameIndex>min(ValidFrameIndices))
        %FrameIndex=FrameIndex-1;
        FrameIndex=ValidFrameIndices(find(ValidFrameIndices < FrameIndex, 1, 'last'));%FrameIndex+1;
        if isempty(BadIndices{FrameIndex})
            BadIndices{FrameIndex} = [];
        end
        
        if isempty(GoodIndices{FrameIndex})
            GoodIndices{FrameIndex} = [];
        end
        
        if isempty(MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = [];
        end
        
        if isempty(BestIndices{FrameIndex})
            BestIndices{FrameIndex} = [];
        end
        
        if isempty(AggregateIndices{FrameIndex})
            AggregateIndices{FrameIndex} = [];
        end
        SpotIndex = 1;
    elseif (ct~=0)&((currentCharacter=='Q'))
        if ~ismember(SpotIndex, BestIndices{FrameIndex})
            BestIndices{FrameIndex} =[BestIndices{FrameIndex} SpotIndex];
        end
        if ~ismember(SpotIndex, GoodIndices{FrameIndex})
            GoodIndices{FrameIndex} =[GoodIndices{FrameIndex} SpotIndex];
        end
        if ismember(SpotIndex, BadIndices{FrameIndex})
            BadIndices{FrameIndex} = BadIndices{FrameIndex}(BadIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = MultiIndices{FrameIndex}(MultiIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, AggregateIndices{FrameIndex})
            AggregateIndices{FrameIndex} = AggregateIndices{FrameIndex}(AggregateIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, BestIndices{FrameIndex})
            set(TempFigure,'Color','c')
        end
        if (SpotIndex<length(Spots(FrameIndex).Fits))
            SpotIndex=SpotIndex+1;
        end
        GoodIndices{FrameIndex} = unique(GoodIndices{FrameIndex});
        BestIndices{FrameIndex} = unique(BestIndices{FrameIndex});
    elseif (ct~=0)&((currentCharacter=='1'))
        if ~ismember(SpotIndex, Cat1Indices{FrameIndex})
            Cat1Indices{FrameIndex} =[Cat1Indices{FrameIndex} SpotIndex];
            Cat1Indices{FrameIndex} = unique(Cat1Indices{FrameIndex});
        end
        if ismember(SpotIndex, Cat2Indices{FrameIndex})
            Cat2Indices{FrameIndex} = Cat2Indices{FrameIndex}(Cat2Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat3Indices{FrameIndex})
            Cat3Indices{FrameIndex} = Cat3Indices{FrameIndex}(Cat3Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat4Indices{FrameIndex})
            Cat4Indices{FrameIndex} = Cat4Indices{FrameIndex}(Cat4Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat5Indices{FrameIndex})
            Cat5Indices{FrameIndex} = Cat5Indices{FrameIndex}(Cat5Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, AggregateIndices{FrameIndex})
            AggregateIndices{FrameIndex} = AggregateIndices{FrameIndex}(AggregateIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = MultiIndices{FrameIndex}(MultiIndices{FrameIndex} ~= SpotIndex);
        end
        if (SpotIndex<length(Spots(FrameIndex).Fits))
            SpotIndex=SpotIndex+1;
        end
        
    elseif (ct~=0)&((currentCharacter=='2'))
        if ~ismember(SpotIndex, Cat2Indices{FrameIndex})
            Cat2Indices{FrameIndex} =[Cat2Indices{FrameIndex} SpotIndex];
            Cat2Indices{FrameIndex} = unique(Cat2Indices{FrameIndex});
        end
        if ismember(SpotIndex, Cat1Indices{FrameIndex})
            Cat1Indices{FrameIndex} = Cat1Indices{FrameIndex}(Cat1Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat3Indices{FrameIndex})
            Cat3Indices{FrameIndex} = Cat3Indices{FrameIndex}(Cat3Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat4Indices{FrameIndex})
            Cat4Indices{FrameIndex} = Cat4Indices{FrameIndex}(Cat4Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat5Indices{FrameIndex})
            Cat5Indices{FrameIndex} = Cat5Indices{FrameIndex}(Cat5Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, AggregateIndices{FrameIndex})
            AggregateIndices{FrameIndex} = AggregateIndices{FrameIndex}(AggregateIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = MultiIndices{FrameIndex}(MultiIndices{FrameIndex} ~= SpotIndex);
        end
        if (SpotIndex<length(Spots(FrameIndex).Fits))
            SpotIndex=SpotIndex+1;
        end
     elseif (ct~=0)&((currentCharacter=='3'))
        if ~ismember(SpotIndex, Cat3Indices{FrameIndex})
            Cat3Indices{FrameIndex} =[Cat3Indices{FrameIndex} SpotIndex];
            Cat3Indices{FrameIndex} = unique(Cat3Indices{FrameIndex});
        end
        if ismember(SpotIndex, Cat1Indices{FrameIndex})
            Cat1Indices{FrameIndex} = Cat1Indices{FrameIndex}(Cat1Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat2Indices{FrameIndex})
            Cat2Indices{FrameIndex} = Cat2Indices{FrameIndex}(Cat2Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat4Indices{FrameIndex})
            Cat4Indices{FrameIndex} = Cat4Indices{FrameIndex}(Cat4Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat5Indices{FrameIndex})
            Cat5Indices{FrameIndex} = Cat5Indices{FrameIndex}(Cat5Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, AggregateIndices{FrameIndex})
            AggregateIndices{FrameIndex} = AggregateIndices{FrameIndex}(AggregateIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = MultiIndices{FrameIndex}(MultiIndices{FrameIndex} ~= SpotIndex);
        end
        if (SpotIndex<length(Spots(FrameIndex).Fits))
            SpotIndex=SpotIndex+1;
        end
    elseif (ct~=0)&((currentCharacter=='4'))
        if ~ismember(SpotIndex, Cat4Indices{FrameIndex})
            Cat4Indices{FrameIndex} =[Cat4Indices{FrameIndex} SpotIndex];
            Cat4Indices{FrameIndex} = unique(Cat4Indices{FrameIndex});
        end
        if ismember(SpotIndex, Cat1Indices{FrameIndex})
            Cat1Indices{FrameIndex} = Cat1Indices{FrameIndex}(Cat1Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat2Indices{FrameIndex})
            Cat2Indices{FrameIndex} = Cat2Indices{FrameIndex}(Cat2Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat3Indices{FrameIndex})
            Cat3Indices{FrameIndex} = Cat3Indices{FrameIndex}(Cat3Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat5Indices{FrameIndex})
            Cat5Indices{FrameIndex} = Cat5Indices{FrameIndex}(Cat5Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, AggregateIndices{FrameIndex})
            AggregateIndices{FrameIndex} = AggregateIndices{FrameIndex}(AggregateIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = MultiIndices{FrameIndex}(MultiIndices{FrameIndex} ~= SpotIndex);
        end
        if (SpotIndex<length(Spots(FrameIndex).Fits))
            SpotIndex=SpotIndex+1;
        end
    elseif (ct~=0)&((currentCharacter=='5'))
        if ~ismember(SpotIndex, Cat5Indices{FrameIndex})
            Cat5Indices{FrameIndex} =[Cat5Indices{FrameIndex} SpotIndex];
            Cat5Indices{FrameIndex} = unique(Cat5Indices{FrameIndex});
        end
        if ismember(SpotIndex, Cat1Indices{FrameIndex})
            Cat1Indices{FrameIndex} = Cat1Indices{FrameIndex}(Cat1Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat2Indices{FrameIndex})
            Cat2Indices{FrameIndex} = Cat2Indices{FrameIndex}(Cat2Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat3Indices{FrameIndex})
            Cat3Indices{FrameIndex} = Cat3Indices{FrameIndex}(Cat3Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat4Indices{FrameIndex})
            Cat4Indices{FrameIndex} = Cat4Indices{FrameIndex}(Cat4Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, AggregateIndices{FrameIndex})
            AggregateIndices{FrameIndex} = AggregateIndices{FrameIndex}(AggregateIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = MultiIndices{FrameIndex}(MultiIndices{FrameIndex} ~= SpotIndex);
        end
        if (SpotIndex<length(Spots(FrameIndex).Fits))
            SpotIndex=SpotIndex+1;
        end
    elseif (ct~=0)&((currentCharacter=='b'))
        CurrentSpot = SpotIndex;
        CurrentFrame = FrameIndex;
        close all
        FigAx = cell(1, MaxStackSize);
        
        TempFigure = figure(1);
        
        set(TempFigure,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
        set(gcf,'color','w');
        SnippetImage2 = zeros(2*snippet_size+1,2*snippet_size+1, MaxStackSize, 'double');
        SnippetImage3 = zeros(2*max_snippet_size+1,2*max_snippet_size+1, MaxStackSize, 'double');
        SnippetImage4 = zeros(2*zoomed_out_snippet_size+1,2*zoomed_out_snippet_size+1, MaxStackSize, 'double');
        clim = [0, MaxPixelValue];
        
        for SubplotIndex = 1:MaxStackSize
            
            if SubplotIndex < MaxStackSize-1
                FigAx{SubplotIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotIndex);
                
                imagesc(FigAx{SubplotIndex} , squeeze(SnippetImage2(:,:,1)), clim);
                colormap(FigAx{SubplotIndex},'gray')
            elseif SubplotIndex == MaxStackSize -1
                colormap(FigAx{SubplotIndex},'gray')%,'parula')
                FigAx{SubplotIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotIndex);
                
                imagesc(FigAx{SubplotIndex} , max(SnippetImage3, [], 3), clim);
            else
                colormap(FigAx{SubplotIndex},'parula')%,'parula')
                FigAx{SubplotIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotIndex);
                
                imagesc(FigAx{SubplotIndex} , max(SnippetImage4, [], 3), clim);
            end
            axis off
        end
        
        for SubplotIndex =1:MaxStackSize
            pos = get(FigAx{SubplotIndex} , 'position');
            pos(1) = SubplotPositionInfo.SubplotXPositions(SubplotIndex);
            pos(3) = SubplotPositionInfo.SubplotWidth;
            pos(2) = SubplotPositionInfo.SubplotYPositions(SubplotIndex);
            pos(4) = SubplotPositionInfo.SubplotHeight;
            set(FigAx{SubplotIndex}, 'position', pos);
        end
        CurrentZ = Spots(FrameIndex).Fits(SpotIndex).z(1);
        plot_title= {strrep(TifList(FrameIndex).name, '_', '\_'), ['Spot Index : ', num2str(SpotIndex)]};
        
        MiddleZ = uint8(floor(median(Spots(FrameIndex).Fits(SpotIndex).z)));
        RefIndex = find(Spots(FrameIndex).Fits(SpotIndex).z == MiddleZ);
        %%
        xmin = Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex)-snippet_size;
        xmax = Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex)+snippet_size;
        ymin = Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex)-snippet_size;
        ymax = Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex)+snippet_size;
        if xmin> 0 & ymin > 0 & xmax <= xDim & ymax <= yDim
            SnippetImage2(:,:, 1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = ...
                MovieCells{FrameIndex}(ymin:ymax,xmin:xmax, Spots(FrameIndex).Fits(SpotIndex).z);
        else
            if xmin <= 0
                xSnip = 1-(xmin);
            else
                xSnip =1;
            end
            if ymin <= 0
                ySnip = 1-(ymin);
            else
                ySnip=1;
            end
            xmaxSnip = xSnip + (xmax-xmin);
            ymaxSnip = ySnip + (ymax-ymin);
            SnippetImage2(ySnip:ymaxSnip, xSnip:xmaxSnip, 1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = MovieCells{FrameIndex}(ymin:ymax,xmin:xmax, Spots(FrameIndex).Fits(SpotIndex).z);
        end
        SnippetImage3 = zeros(2*max_snippet_size+1,2*max_snippet_size+1,MaxStackSize, 'double');
        SnippetImage4 = zeros(2*zoomed_out_snippet_size+1,2*zoomed_out_snippet_size+1,MaxStackSize, 'double');
        
        max_xmin = double(Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex))-max_snippet_size;
        max_xmax = double(Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex))+max_snippet_size;
        max_ymin = double(Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex))-max_snippet_size;
        max_ymax = double(Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex))+max_snippet_size;
        if max_xmin> 0 & max_ymin > 0 & max_xmax <= xDim & max_ymax <= yDim
            SnippetImage3 = MovieCells{FrameIndex}(max_ymin:max_ymax,max_xmin:max_xmax,Spots(FrameIndex).Fits(SpotIndex).z);%.*double(SpotMat(xmin:xmax,ymin:ymax,ReorderedIndex));
        else
            if max_xmin <= 0
                max_xSnip = 1-(max_xmin);
                max_xmin = 1;
            else
                max_xSnip =1;
            end
            if max_ymin <= 0
                max_ySnip = 1-(max_ymin);
                max_ymin = 1;
            else
                max_ySnip=1;
            end
            if max_xmax > xDim
                max_xmax = xDim;
            end
            if max_ymax > yDim
                max_ymax = yDim;
            end
            max_xmaxSnip = max_xSnip + (max_xmax-max_xmin);
            max_ymaxSnip = max_ySnip + (max_ymax-max_ymin);
            SnippetImage3(max_ySnip:max_ymaxSnip,max_xSnip:max_xmaxSnip,  1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = MovieCells{FrameIndex}(max_ymin:max_ymax,max_xmin:max_xmax,Spots(FrameIndex).Fits(SpotIndex).z);%.*double(SpotMat(xmin:xmax,ymin:ymax,ReorderedIndex));
        end
        
        
        zoomed_xmin = double(Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex))-zoomed_out_snippet_size;
        zoomed_xmax = double(Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex))+zoomed_out_snippet_size;
        zoomed_ymin = double(Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex))-zoomed_out_snippet_size;
        zoomed_ymax = double(Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex))+zoomed_out_snippet_size;
        if zoomed_xmin> 0 & zoomed_ymin > 0 & zoomed_xmax <= xDim & zoomed_ymax <= yDim
            SnippetImage4(:,:, 1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = ...
                MovieCells{FrameIndex}(zoomed_ymin:zoomed_ymax,zoomed_xmin:zoomed_xmax, Spots(FrameIndex).Fits(SpotIndex).z);
        else
            if zoomed_xmin <= 0
                zoomed_xSnip = 1-(zoomed_xmin);
                zoomed_xmin = 1;
            else
                zoomed_xSnip =1;
            end
            if zoomed_ymin <= 0
                zoomed_ySnip = 1-(zoomed_ymin);
                zoomed_ymin = 1;
            else
                zoomed_ySnip=1;
            end
            
            
            if zoomed_xmax > xDim
                zoomed_xmax = xDim;
            end
            if zoomed_ymax > yDim
                zoomed_ymax = yDim;
            end
            
            zoomed_xmaxSnip = zoomed_xSnip + (zoomed_xmax-zoomed_xmin);
            zoomed_ymaxSnip = zoomed_ySnip + (zoomed_ymax-zoomed_ymin);
            SnippetImage4(zoomed_ySnip:zoomed_ymaxSnip, zoomed_xSnip:zoomed_xmaxSnip, 1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = MovieCells{FrameIndex}(zoomed_ymin:zoomed_ymax,zoomed_xmin:zoomed_xmax, Spots(FrameIndex).Fits(SpotIndex).z);
        end
        
        
        ScatterHandlesRed =cell(1,MaxStackSize);
        ScatterHandlesYellow =cell(1,MaxStackSize);
        ScatterHandlesBlue =cell(1,MaxStackSize);
        for SubplotIndex = 1:MaxStackSize-2
            set(TempFigure, 'currentaxes', FigAx{SubplotIndex});
            set( FigAx{SubplotIndex}, 'YLimMode', 'manual' )
            set( FigAx{SubplotIndex}, 'XLimMode', 'manual' )
            xlim([0.5 2*snippet_size+1.5]);
            ylim([0.5 2*snippet_size+1.5]);
            if SubplotIndex <= length(Spots(FrameIndex).Fits(SpotIndex).z)
                
                imagesc(FigAx{SubplotIndex} , squeeze(SnippetImage2(:,:,SubplotIndex)), clim);
                title(['z = ', num2str(Spots(FrameIndex).Fits(SpotIndex).z(SubplotIndex))])
                hold on
                colormap(FigAx{SubplotIndex},'gray')
                axis off
                
                plot(FigAx{SubplotIndex}, pline_x, pline_y, '.r');
                row_pixels= Spots(FrameIndex).Fits(SpotIndex).PixelList{SubplotIndex}(:,1).'-double(xmin)+1;
                col_pixels= Spots(FrameIndex).Fits(SpotIndex).PixelList{SubplotIndex}(:,2).'-double(ymin)+1;
                ScatterHandlesRed{SubplotIndex} = scatter(row_pixels,col_pixels, 'r*', 'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
                [mid_row, mid_cols] = find(squeeze(SnippetImage2(:,:,SubplotIndex)) >= MaxPixelValue/2);
                ScatterHandlesYellow{SubplotIndex} = scatter(mid_cols,mid_row,  'y*',  'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
                [high_row, high_cols] = find(squeeze(SnippetImage2(:,:,SubplotIndex)) >= MaxPixelValue);
                ScatterHandlesBlue{SubplotIndex} = scatter(high_cols,high_row,  'c*',  'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
                
                hold off
            else
                imagesc(FigAx{SubplotIndex}, zeros(2*snippet_size+1,2*snippet_size+1,'double'), clim);
                hold on
                plot(FigAx{SubplotIndex}, pline_x, pline_y, '.r');
                [mid_row, mid_cols] = find(squeeze(SnippetImage2(:,:,SubplotIndex)) >= MaxPixelValue/2);
                [high_row, high_cols] = find(squeeze(SnippetImage2(:,:,SubplotIndex)) >= MaxPixelValue);
                ScatterHandlesRed{SubplotIndex} = scatter(high_cols,high_row, 'r*', 'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
                
                ScatterHandlesYellow{SubplotIndex} = scatter(mid_cols, mid_row,  'y*',  'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
                ScatterHandlesBlue{SubplotIndex} = scatter(high_cols, high_row,  'c*',  'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
                colormap(FigAx{SubplotIndex},'gray')
                axis off
                hold off
                title('')
            end
            
            
            
        end
        
        set(TempFigure, 'currentaxes', FigAx{MaxStackSize-1});
        set( FigAx{MaxStackSize-1}, 'YLimMode', 'manual' )
        set( FigAx{MaxStackSize-1}, 'XLimMode', 'manual' )
        xlim([0.5 2*max_snippet_size+1.5]);
        ylim([0.5 2*max_snippet_size+1.5]);
        imagesc(FigAx{MaxStackSize-1} , sum(SnippetImage3,3), clim*1);% max(SnippetImage3,[],3), clim);
        hold on
        
        title('Max Proj.')
        colormap(FigAx{MaxStackSize-1},'gray')%,'parula')
        axis off
        plot(FigAx{MaxStackSize-1}, max_pline_x, max_pline_y, '.r');
        hold off
        
        set(TempFigure, 'currentaxes', FigAx{MaxStackSize});
        set( FigAx{MaxStackSize}, 'YLimMode', 'manual' )
        set( FigAx{MaxStackSize}, 'XLimMode', 'manual' )
        xlim([0.5 2*zoomed_out_snippet_size+1.5]);
        ylim([0.5 2*zoomed_out_snippet_size+1.5]);
        imagesc(FigAx{MaxStackSize} , sum(SnippetImage4,3), clim*1);% max(SnippetImage3,[],3), clim);
        hold on
        
        title('Max Proj.')
        colormap(FigAx{MaxStackSize},'parula')%,'parula')
        axis off
        plot(FigAx{MaxStackSize}, zoomed_pline_x, zoomed_pline_y,'.r',...
            'MarkerSize',5);
        %         rowVector = Spots(FrameIndex).Fits(SpotIndex).xDoG-double(zoomed_xmin)+1;
        %         colVector = Spots(FrameIndex).Fits(SpotIndex).yDoG-double(zoomed_ymin)+1;
        %
        %
        %         ScatterHandlesRed{MaxStackSize} = scatter(FigAx{MaxStackSize}, rowVector, colVector,10, 'MarkerFaceColor','r','MarkerEdgeColor','r');
        hold off
        
        
        sgtitle(plot_title)
        
        if isempty(BadIndices{FrameIndex})
            BadIndices{FrameIndex} = [];
        end
        
        if isempty(GoodIndices{FrameIndex})
            GoodIndices{FrameIndex} = [];
        end
        
        if isempty(MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = [];
        end
        
        if isempty(BestIndices{FrameIndex})
            BestIndices{FrameIndex} = [];
        end
        
        if isempty(AggregateIndices{FrameIndex})
            AggregateIndices{FrameIndex} = [];
        end
        if isempty(Cat1Indices{FrameIndex})
            Cat1Indices{FrameIndex} = [];
        end
        if isempty(Cat2Indices{FrameIndex})
            Cat2Indices{FrameIndex} = [];
        end
        if isempty(Cat3Indices{FrameIndex})
            Cat3Indices{FrameIndex} = [];
        end
        if isempty(Cat4Indices{FrameIndex})
            Cat4Indices{FrameIndex} = [];
        end
        if isempty(Cat5Indices{FrameIndex})
            Cat5Indices{FrameIndex} = [];
        end
        
        
        
        if ismember(SpotIndex, BestIndices{FrameIndex})
            set(TempFigure,'Color','c')
        elseif ismember(SpotIndex, AggregateIndices{FrameIndex})
            set(TempFigure,'Color','m')
        elseif ismember(SpotIndex, BadIndices{FrameIndex})
            set(TempFigure,'Color','r')
        elseif ismember(SpotIndex, GoodIndices{FrameIndex})
            set(TempFigure,'Color','g')
        elseif ismember(SpotIndex, MultiIndices{FrameIndex})
            set(TempFigure,'Color','b')
        else
            set(TempFigure,'Color','y')
        end
        
        hold off
        
        %%
        SpotIndex=CurrentSpot;
        FrameIndex= CurrentFrame ;
        currentCharacter = 0;
        set(0, 'CurrentFigure', TempFigure)
        
    elseif (ct~=0)&((currentCharacter=='y') | (currentCharacter=='q'))
        GoodIndices{FrameIndex} =[GoodIndices{FrameIndex} SpotIndex];
        if ismember(SpotIndex, BadIndices{FrameIndex})
            BadIndices{FrameIndex} = BadIndices{FrameIndex}(BadIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = MultiIndices{FrameIndex}(MultiIndices{FrameIndex} ~= SpotIndex);
        end
        if (SpotIndex<length(Spots(FrameIndex).Fits))
            SpotIndex=SpotIndex+1;
        end
        if ismember(SpotIndex, GoodIndices{FrameIndex})
            set(TempFigure,'Color','g')
        end
    elseif (ct~=0)&((currentCharacter=='n') | (currentCharacter=='w'))
        BadIndices{FrameIndex} =[BadIndices{FrameIndex} SpotIndex];
        if ismember(SpotIndex, GoodIndices{FrameIndex})
            GoodIndices{FrameIndex} = GoodIndices{FrameIndex}(GoodIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = MultiIndices{FrameIndex}(MultiIndices{FrameIndex} ~= SpotIndex);
        end
        if (SpotIndex<length(Spots(FrameIndex).Fits))
            SpotIndex=SpotIndex+1;
        end
        if ismember(SpotIndex, BadIndices{FrameIndex})
            set(TempFigure,'Color','r')
        end
    elseif (ct~=0)&((currentCharacter=='W'))
        AggregateIndices{FrameIndex} =[AggregateIndices{FrameIndex} SpotIndex];
        AggregateIndices{FrameIndex} = unique(AggregateIndices{FrameIndex});
        BadIndices{FrameIndex} =[BadIndices{FrameIndex} SpotIndex];
        BadIndices{FrameIndex} = unique(BadIndices{FrameIndex});
        if ismember(SpotIndex, GoodIndices{FrameIndex})
            GoodIndices{FrameIndex} = GoodIndices{FrameIndex}(GoodIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, BestIndices{FrameIndex})
            BestIndices{FrameIndex} = BestIndices{FrameIndex}(BestIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = MultiIndices{FrameIndex}(MultiIndices{FrameIndex} ~= SpotIndex);
        end
        if (SpotIndex<length(Spots(FrameIndex).Fits))
            SpotIndex=SpotIndex+1;
        end
        if ismember(SpotIndex, AggregateIndices{FrameIndex})
            set(TempFigure,'Color','m')
        end
    elseif (ct~=0)&(currentCharacter=='j')
        try
            iJump = inputdlg('Jump to Spot Index:', ...
                'Spot Index:');
            iJump = str2double(iJump{1});
        catch
            iJump = 0;
        end
        if ( iJump <= length(Spots(FrameIndex).Fits)) & (iJump > 0)
            SpotIndex = iJump;
        elseif ( iJump <= length(Spots(FrameIndex).Fits))
            SpotIndex = 1;
        else
            SpotIndex = length(Spots(FrameIndex).Fits);
        end
    elseif (ct~=0)&((currentCharacter=='o')  | (currentCharacter=='e'))
        MultiIndices{FrameIndex} =[MultiIndices{FrameIndex} SpotIndex];
        MultiIndices{FrameIndex} = unique(MultiIndices{FrameIndex});
        if ismember(SpotIndex, GoodIndices{FrameIndex})
            GoodIndices{FrameIndex} = GoodIndices{FrameIndex}(GoodIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, BadIndices{FrameIndex})
            BadIndices{FrameIndex} = BadIndices{FrameIndex}(BadIndices{FrameIndex} ~= SpotIndex);
        end
        if (SpotIndex<length(Spots(FrameIndex).Fits))
            SpotIndex=SpotIndex+1;
        end
        if ismember(SpotIndex, MultiIndices{FrameIndex})
            set(TempFigure,'Color','b')
        end
    elseif (ct~=0)&(currentCharacter=='r')
        if ismember(SpotIndex, MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = MultiIndices{FrameIndex}(MultiIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, BestIndices{FrameIndex})
            BestIndices{FrameIndex} = BestIndices{FrameIndex}(BestIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, AggregateIndices{FrameIndex})
            AggregateIndices{FrameIndex} = AggregateIndices{FrameIndex}(AggregateIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, GoodIndices{FrameIndex})
            GoodIndices{FrameIndex} = GoodIndices{FrameIndex}(GoodIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, BadIndices{FrameIndex})
            BadIndices{FrameIndex} = BadIndices{FrameIndex}(BadIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat1Indices{FrameIndex})
            Cat1Indices{FrameIndex} = Cat1Indices{FrameIndex}(Cat1Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat2Indices{FrameIndex})
            Cat2Indices{FrameIndex} = Cat2Indices{FrameIndex}(Cat2Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat3Indices{FrameIndex})
            Cat3Indices{FrameIndex} = Cat3Indices{FrameIndex}(Cat3Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat4Indices{FrameIndex})
            Cat4Indices{FrameIndex} = Cat4Indices{FrameIndex}(Cat4Indices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, Cat5Indices{FrameIndex})
            Cat5Indices{FrameIndex} = Cat5Indices{FrameIndex}(Cat5Indices{FrameIndex} ~= SpotIndex);
        end
    elseif (ct~=0)&(currentCharacter=='p')
        for fr = 1:NumFrames
            MaxLabeledIndex = max([BadIndices{fr} GoodIndices{fr} AggregateIndices{fr} BestIndices{fr} MultiIndices{fr}]);
            disp(['Number of good spots in frame ',num2str(fr),': ',num2str(length(GoodIndices{fr} ))])
            disp(['Max labeled spot in frame ',num2str(fr),': ',num2str(MaxLabeledIndex), ' (out of ',num2str(length(Spots(fr).Fits)),')']);
        end
    elseif (ct~=0)&(currentCharacter=='s')
        GoodSpots = Spots;
        BadSpots= Spots;
        MultiSpots = Spots;
        BestSpots = Spots;
        AggregateSpots = Spots;
        Cat1Spots = Spots;
        Cat2Spots = Spots;
        Cat3Spots = Spots;
        Cat4Spots = Spots;
        Cat5Spots = Spots;
        for i = 1:length(Spots)
            GoodIndices{i} = unique(GoodIndices{i});
            BadIndices{i} = unique(BadIndices{i});
            MultiIndices{i} = unique(MultiIndices{i});
            BestIndices{i} = unique(BestIndices{i});
            AggregateIndices{i} = unique(AggregateIndices{i});
            Cat1Indices{i}= unique(Cat1Indices{i});
            Cat2Indices{i}= unique(Cat2Indices{i});
            Cat3Indices{i}= unique(Cat3Indices{i});
            Cat4Indices{i}= unique(Cat4Indices{i});
            Cat5Indices{i}= unique(Cat5Indices{i});
            GoodSpots(i).Fits = Spots(i).Fits(GoodIndices{i}) ;
            BestSpots(i).Fits = Spots(i).Fits(BestIndices{i}) ;
            BadSpots(i).Fits = Spots(i).Fits(BadIndices{i}) ;
            MultiSpots(i).Fits = Spots(i).Fits(MultiIndices{i}) ;
            AggregateSpots(i).Fits = Spots(i).Fits(AggregateIndices{i}) ;
            Cat1Spots(i).Fits = Spots(i).Fits(Cat1Indices{i});
            Cat2Spots(i).Fits = Spots(i).Fits(Cat2Indices{i});
            Cat3Spots(i).Fits = Spots(i).Fits(Cat3Indices{i});
            Cat4Spots(i).Fits = Spots(i).Fits(Cat4Indices{i});
            Cat5Spots(i).Fits = Spots(i).Fits(Cat5Indices{i});
        end
        
        
        outpath = [liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat'];
        save(outpath, 'GoodSpots', 'BadSpots', 'MultiSpots', 'BestSpots', 'AggregateSpots',...
            'GoodIndices', 'BadIndices', 'MultiIndices', 'BestIndices', 'AggregateIndices',...
            'Cat1Spots','Cat2Spots','Cat3Spots','Cat4Spots','Cat5Spots',...
            'Cat1Indices','Cat2Indices','Cat3Indices','Cat4Indices','Cat5Indices');
        
        
    end
end
%%
close all

for fr = 1:NumFrames
    MaxLabeledIndex = max([BadIndices{fr} GoodIndices{fr} AggregateIndices{fr} BestIndices{fr} MultiIndices{fr}]);
    disp(['Number of good spots in frame ',num2str(fr),': ',num2str(length(GoodIndices{fr} ))])
    disp(['Max labeled spot in frame ',num2str(fr),': ',num2str(MaxLabeledIndex), ' (out of ',num2str(length(Spots(fr).Fits)),')']);
end

GoodSpots = Spots;
BadSpots= Spots;
MultiSpots = Spots;
BestSpots = Spots;
AggregateSpots = Spots;
Cat1Spots = Spots;
        Cat2Spots = Spots;
        Cat3Spots = Spots;
        Cat4Spots = Spots;
        Cat5Spots = Spots;
for i = 1:length(Spots)
    GoodIndices{i} = unique(GoodIndices{i});
    BadIndices{i} = unique(BadIndices{i});
    MultiIndices{i} = unique(MultiIndices{i});
    BestIndices{i} = unique(BestIndices{i});
    AggregateIndices{i} = unique(AggregateIndices{i});
    GoodSpots(i).Fits = Spots(i).Fits(GoodIndices{i}) ;
    BadSpots(i).Fits = Spots(i).Fits(BadIndices{i}) ;
    MultiSpots(i).Fits = Spots(i).Fits(MultiIndices{i}) ;
    BestSpots(i).Fits = Spots(i).Fits(BestIndices{i}) ;
    AggregateSpots(i).Fits = Spots(i).Fits(AggregateIndices{i}) ;
    Cat1Spots(i).Fits = Spots(i).Fits(Cat1Indices{i});
            Cat2Spots(i).Fits = Spots(i).Fits(Cat2Indices{i});
            Cat3Spots(i).Fits = Spots(i).Fits(Cat3Indices{i});
            Cat4Spots(i).Fits = Spots(i).Fits(Cat4Indices{i});
            Cat5Spots(i).Fits = Spots(i).Fits(Cat5Indices{i});
end


outpath = [liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat'];
save(outpath, 'GoodSpots', 'BadSpots', 'MultiSpots', 'BestSpots', 'AggregateSpots',...
    'GoodIndices', 'BadIndices', 'MultiIndices', 'BestIndices', 'AggregateIndices',...
            'Cat1Spots','Cat2Spots','Cat3Spots','Cat4Spots','Cat5Spots',...
            'Cat1Indices','Cat2Indices','Cat3Indices','Cat4Indices','Cat5Indices');

NumGoodSpots = length(GoodIndices{FrameIndex});
%%