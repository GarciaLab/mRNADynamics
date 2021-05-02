function [FrameIndex, SpotIndex,NumFrames, NumGoodSpots] = Store3DNanocageInfo(Prefix, MaxPixelValue, FrameIndex, SpotIndex)
if ~exist('MaxPixelValue', 'var')
    MaxPixelValue = 20;
end
liveExperiment = LiveExperiment(Prefix);
Spots = getSpots(liveExperiment);
snippet_size = ceil(1300/(2*liveExperiment.pixelSize_nm));
xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
NumFrames = length(Spots);
TifList = dir([liveExperiment.preFolder, filesep, '*ch01.tif']);
ProbList = dir([liveExperiment.procFolder,'dogs',  filesep, 'prob*']);
MovieCells = cell(1, NumFrames);
ProbCells = cell(1, NumFrames);
NumSlices = 0;
for i = 1:NumFrames
    MovieCells{i} = imreadStack([TifList(i).folder, filesep,  TifList(i).name]);
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
else
    GoodIndices = cell(1, length(Spots));
    BadIndices = cell(1, length(Spots));
    MultiIndices = cell(1, length(Spots));
    BestIndices = cell(1, length(Spots));
    AggregateIndices = cell(1, length(Spots));
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
theta = 0 : (2 * pi / 10000) : (2 * pi);
pline_x = snippet_size/2 * cos(theta) + snippet_size+1;
pline_y = snippet_size/2  * sin(theta) + snippet_size+1;

close all
FigAx = cell(1, MaxStackSize);

TempFigure = figure(1);

set(TempFigure,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
set(gcf,'color','w');
SnippetImage2 = zeros(2*snippet_size+1,2*snippet_size+1, MaxStackSize, 'double');

clim = [0, MaxPixelValue];

for SubplotIndex = 1:MaxStackSize
    FigAx{SubplotIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotIndex);
    
    imagesc(FigAx{SubplotIndex} , squeeze(SnippetImage2(:,:,1)), clim);
    colormap('gray')
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

ScatterHandlesRed =cell(1,MaxStackSize);
ScatterHandlesYellow =cell(1,MaxStackSize);
ScatterHandlesBlue =cell(1,MaxStackSize);
for SubplotIndex = 1:MaxStackSize
    set(TempFigure, 'currentaxes', FigAx{SubplotIndex});
    set( FigAx{SubplotIndex}, 'YLimMode', 'manual' )
    set( FigAx{SubplotIndex}, 'XLimMode', 'manual' )
    xlim([0.5 2*snippet_size+1.5]);
    ylim([0.5 2*snippet_size+1.5]);
    if SubplotIndex <= length(Spots(FrameIndex).Fits(SpotIndex).z)
        
        imagesc(FigAx{SubplotIndex} , squeeze(SnippetImage2(:,:,SubplotIndex)), clim);
        title(['z = ', num2str(Spots(FrameIndex).Fits(SpotIndex).z(SubplotIndex))])
        hold on
        colormap('gray')
        axis off
        
        plot(FigAx{SubplotIndex}, pline_x, pline_y, '.r');
        row_pixels= Spots(FrameIndex).Fits(SpotIndex).PixelList{SubplotIndex}(:,1).'-double(xmin)+1;
        col_pixels= Spots(FrameIndex).Fits(SpotIndex).PixelList{SubplotIndex}(:,2).'-double(ymin)+1;
        ScatterHandlesRed{SubplotIndex} = scatter(row_pixels,col_pixels, 'r*', 'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
        [high_row, high_cols] = find(squeeze(SnippetImage2(:,:,SubplotIndex)) > MaxPixelValue);
        ScatterHandlesYellow{SubplotIndex} = scatter(high_cols,high_row,  'y*',  'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
        
        hold off
    else
        imagesc(FigAx{SubplotIndex}, zeros(2*snippet_size+1,2*snippet_size+1,'double'), clim);
        hold on
        plot(FigAx{SubplotIndex}, pline_x, pline_y, '.r');
        [high_row, high_cols] = find(squeeze(SnippetImage2(:,:,SubplotIndex)) > MaxPixelValue);
        ScatterHandlesRed{SubplotIndex} = scatter(high_cols,high_row, 'r*', 'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
        
        ScatterHandlesYellow{SubplotIndex} = scatter(high_cols, high_row,  'y*',  'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
        colormap('gray')
        axis off
        hold off
        title('')
    end
    
    
    
end


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
    plot_title= {[strrep(TifList(FrameIndex).name, '_', '\_'), ' (Frame Count: ', num2str(length(Spots)), ')' ], ['Spot Index: ',...
        num2str(SpotIndex),'/',num2str(length(Spots(FrameIndex).Fits)),...
        ', Z Count: ', num2str(length(Spots(FrameIndex).Fits(SpotIndex).z)),...
        ', Good Spots: ', num2str(NumGoodSpots), ', Bad Spots: ',...
        num2str(NumBadSpots), ', Multi Spots: ', num2str(NumMultiSpots),...
        ', Best Spots: ', num2str(NumBestSpots),...
        ', Aggregate Spots: ', num2str(NumAggregates)]};
    
    MiddleZ = uint8(floor(median(Spots(FrameIndex).Fits(SpotIndex).z)));
    RefIndex = find(Spots(FrameIndex).Fits(SpotIndex).z == MiddleZ);
    %Load subsequent images
    
    
    
    
    SnippetImage = zeros(2*snippet_size+1,2*snippet_size+1,MaxStackSize, 'double');
    xmin = Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex)-snippet_size;
    xmax = Spots(FrameIndex).Fits(SpotIndex).xDoG(RefIndex)+snippet_size;
    ymin = Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex)-snippet_size;
    ymax = Spots(FrameIndex).Fits(SpotIndex).yDoG(RefIndex)+snippet_size;
    if xmin> 0 & ymin > 0 & xmax <= xDim & ymax <= yDim
        SnippetImage = MovieCells{FrameIndex}(ymin:ymax,xmin:xmax,Spots(FrameIndex).Fits(SpotIndex).z);%.*double(SpotMat(xmin:xmax,ymin:ymax,ReorderedIndex));
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
        SnippetImage(ySnip:ymaxSnip,xSnip:xmaxSnip,  1:length(Spots(FrameIndex).Fits(SpotIndex).z)) = MovieCells{FrameIndex}(ymin:ymax,xmin:xmax,Spots(FrameIndex).Fits(SpotIndex).z);%.*double(SpotMat(xmin:xmax,ymin:ymax,ReorderedIndex));
    end
    
    for SubplotIndex = 1:MaxStackSize
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
            [high_row, high_cols] = find(squeeze(SnippetImage(:,:,SubplotIndex)) > MaxPixelValue);
            set(ScatterHandlesYellow{SubplotIndex},'XData', high_cols,'YData', high_row);
            set(FigAx{SubplotIndex}.Children(end-3),'Visible','on'); %'off' or 'on'

        else
            FigAx{SubplotIndex}.Children(end).CData = zeros(2*snippet_size+1,2*snippet_size+1,'double');
            FigAx{SubplotIndex}.Title.String = '';
            
            %[high_row, high_cols] = find(FigAx{SubplotIndex}.Children(end).CData  > MaxPixelValue);
            %set(ScatterHandlesRed{SubplotIndex},'XData', high_cols,'YData', high_row);
            %set(ScatterHandlesYellow{SubplotIndex},'XData', high_cols,'YData', high_row);
            set(FigAx{SubplotIndex}.Children(end-2),'Visible','off'); %'off' or 'on'
            set(FigAx{SubplotIndex}.Children(end-3),'Visible','off'); %'off' or 'on'
        end
        %hold(FigAx{SubplotIndex}, 'off')
        
        
    end
    
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
        
    elseif (ct~=0)&(currentCharacter=='>')&(FrameIndex<max(ValidFrameIndices))
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
        
    elseif (ct~=0)&(currentCharacter=='<')&(FrameIndex>min(ValidFrameIndices))
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
        BestIndices{FrameIndex} =[BestIndices{FrameIndex} SpotIndex];
        GoodIndices{FrameIndex} =[GoodIndices{FrameIndex} SpotIndex];
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
        
    elseif (ct~=0)&((currentCharacter=='b'))
        CurrentSpot = SpotIndex;
        CurrentFrame = FrameIndex;
        close all
        FigAx = cell(1, MaxStackSize);
        
        TempFigure = figure(1);
        
        set(TempFigure,'units', 'normalized', 'position',[0.05, 0.05, SubplotPositionInfo.SubFigDims(1), SubplotPositionInfo.SubFigDims(2)]);
        set(gcf,'color','w');
        SnippetImage2 = zeros(2*snippet_size+1,2*snippet_size+1, MaxStackSize, 'double');
        
        clim = [0, MaxPixelValue];
        
        for SubplotIndex = 1:MaxStackSize
            FigAx{SubplotIndex} = subplot(SubplotPositionInfo.SubplotDims(1), SubplotPositionInfo.SubplotDims(2), SubplotIndex);
            
            imagesc(FigAx{SubplotIndex} , squeeze(SnippetImage2(:,:,1)), clim);
            colormap('gray')
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
        
        ScatterHandlesRed =cell(1,MaxStackSize);
        ScatterHandlesYellow =cell(1,MaxStackSize);
        ScatterHandlesBlue =cell(1,MaxStackSize);
        for SubplotIndex = 1:MaxStackSize
            set(TempFigure, 'currentaxes', FigAx{SubplotIndex});
            set( FigAx{SubplotIndex}, 'YLimMode', 'manual' )
            set( FigAx{SubplotIndex}, 'XLimMode', 'manual' )
            xlim([0.5 2*snippet_size+1.5]);
            ylim([0.5 2*snippet_size+1.5]);
            if SubplotIndex <= length(Spots(FrameIndex).Fits(SpotIndex).z)
                
                imagesc(FigAx{SubplotIndex} , squeeze(SnippetImage2(:,:,SubplotIndex)), clim);
                title(['z = ', num2str(Spots(FrameIndex).Fits(SpotIndex).z(SubplotIndex))])
                hold on
                colormap('gray')
                axis off
                
                plot(FigAx{SubplotIndex}, pline_x, pline_y, '.r');
                row_pixels= Spots(FrameIndex).Fits(SpotIndex).PixelList{SubplotIndex}(:,1).'-double(xmin)+1;
                col_pixels= Spots(FrameIndex).Fits(SpotIndex).PixelList{SubplotIndex}(:,2).'-double(ymin)+1;
                ScatterHandlesRed{SubplotIndex} = scatter(row_pixels,col_pixels, 'r*', 'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
                [high_row, high_cols] = find(squeeze(SnippetImage2(:,:,SubplotIndex)) > MaxPixelValue);
                ScatterHandlesYellow{SubplotIndex} = scatter(high_cols,high_row,  'y*',  'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
                
                hold off
            else
                imagesc(FigAx{SubplotIndex}, zeros(2*snippet_size+1,2*snippet_size+1,'double'), clim);
                hold on
                plot(FigAx{SubplotIndex}, pline_x, pline_y, '.r');
                [high_row, high_cols] = find(squeeze(SnippetImage2(:,:,SubplotIndex)) > MaxPixelValue);
                ScatterHandlesRed{SubplotIndex} = scatter(high_cols,high_row, 'r*', 'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
                
                ScatterHandlesYellow{SubplotIndex} = scatter(high_cols, high_row,  'y*',  'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
                colormap('gray')
                axis off
                hold off
                title('')
            end
            
            
            
        end
        
        
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
        BadIndices{FrameIndex} =[BadIndices{FrameIndex} SpotIndex];
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
        for i = 1:length(Spots)
            GoodIndices{i} = unique(GoodIndices{i});
            BadIndices{i} = unique(BadIndices{i});
            MultiIndices{i} = unique(MultiIndices{i});
            BestIndices{i} = unique(BestIndices{i});
            AggregateIndices{i} = unique(AggregateIndices{i});
            GoodSpots(i).Fits = Spots(i).Fits(GoodIndices{i}) ;
            BestSpots(i).Fits = Spots(i).Fits(BestIndices{i}) ;
            BadSpots(i).Fits = Spots(i).Fits(BadIndices{i}) ;
            MultiSpots(i).Fits = Spots(i).Fits(MultiIndices{i}) ;
            AggregateSpots(i).Fits = Spots(i).Fits(AggregateIndices{i}) ;
        end
        
        
        outpath = [liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat'];
        save(outpath, 'GoodSpots', 'BadSpots', 'MultiSpots', 'BestSpots', 'AggregateSpots',...
            'GoodIndices', 'BadIndices', 'MultiIndices', 'BestIndices', 'AggregateIndices');
        
        
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
end


outpath = [liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat'];
save(outpath, 'GoodSpots', 'BadSpots', 'MultiSpots', 'BestSpots', 'AggregateSpots',...
    'GoodIndices', 'BadIndices', 'MultiIndices', 'BestIndices', 'AggregateIndices');

NumGoodSpots = length(GoodIndices{FrameIndex});
%%