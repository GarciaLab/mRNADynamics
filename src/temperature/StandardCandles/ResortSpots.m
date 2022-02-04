
liveExperiment = LiveExperiment(Prefix);
Spots = getSpots(liveExperiment);
snippet_size = ceil(1300/(2*liveExperiment.pixelSize_nm));
max_snippet_size = snippet_size*2;
zoomed_out_snippet_size = snippet_size*6;
xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
NumFrames = length(Spots);
UnsortedSpots = Spots;
if whos(var2str(Spots)).bytes < 2E9
    save([liveExperiment.resultsFolder,...
        filesep, 'UnsortedSpots.mat'], 'UnsortedSpots', '-v6');
else
    save([liveExperiment.resultsFolder,...
        filesep, 'UnsortedSpots.mat'], 'UnsortedSpots', '-v7.3', '-nocompression');
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


SortedGoodIndices = cell(1, length(Spots));
SortedBadIndices = cell(1,length(Spots));
SortedBestIndices = cell(1,length(Spots));
SortedAggregateIndices = cell(1, length(Spots));
SortedMultiIndices = cell(1,length(Spots));


UnsortedGoodIndices = GoodIndices;
UnsortedBadIndices = BadIndices;
UnsortedBestIndices = BestIndices;
UnsortedAggregateIndices = AggregateIndices;
UnsortedMultiIndices = MultiIndices;

%%
FrameIndex = 1;
for FrameIndex = 1:1%length(Spots)

if ~isempty(Spots(FrameIndex).Fits)
    continue
end

NumSpots = length(Spots(FrameIndex).Fits);

Spots(FrameIndex).Fits(1).middleX = NaN;
Spots(FrameIndex).Fits(1).middleY = NaN;
middlex = NaN(1,NumSpots);
middley = NaN(1,NumSpots);
for spotIndex = 1:NumSpots
    MiddleZ = uint8(floor(median(Spots(FrameIndex).Fits(spotIndex).z)));
    RefIndex = find(Spots(FrameIndex).Fits(spotIndex).z == MiddleZ);
    middlex(spotIndex) = Spots(FrameIndex).Fits(spotIndex).xDoG(RefIndex);
    middley(spotIndex) = Spots(FrameIndex).Fits(spotIndex).yDoG(RefIndex);
end
   
midxcell = num2cell(middlex);
midycell = num2cell(middley);
[Spots(FrameIndex).Fits.middleX] = midxcell{:};
[Spots(FrameIndex).Fits.middleY] = midycell{:};

Spots(FrameIndex).Fits(1).xTiles = NaN;
xtcell = num2cell(floor([Spots(FrameIndex).Fits(:).middleX]/zoomed_out_snippet_size)+1);
[Spots(FrameIndex).Fits.xTiles] = xtcell{:};

Spots(FrameIndex).Fits(1).yTiles = NaN;
ytcell = num2cell(floor([Spots(FrameIndex).Fits(:).middleY]/zoomed_out_snippet_size)+1);
[Spots(FrameIndex).Fits.yTiles] = ytcell{:};

PositionalInfo = NaN(NumSpots,5);
PositionalInfo(:,1) = [Spots(FrameIndex).Fits(:).yTiles];
PositionalInfo(:,2) = [Spots(FrameIndex).Fits(:).xTiles];
PositionalInfo(:,3) = [Spots(FrameIndex).Fits(:).middleY];
PositionalInfo(:,4) = [Spots(FrameIndex).Fits(:).middleX];
PositionalInfo(:,5) = [Spots(FrameIndex).Fits(:).FirstZ];

[SortedPositionInfo, SortIndex ] =  sortrows(PositionalInfo, [1,2,5,3,4]);
%%
SortedGoodIndices{FrameIndex} = NaN(1,length(GoodIndices{FrameIndex}));
for spotIndex = 1:length(GoodIndices{FrameIndex})
    SortedGoodIndices{FrameIndex}(spotIndex) = find(SortIndex == GoodIndices{FrameIndex}(spotIndex));
end

SortedBadIndices{FrameIndex} = NaN(1,length(BadIndices{FrameIndex}));
for spotIndex = 1:length(BadIndices{FrameIndex})
    SortedBadIndices{FrameIndex}(spotIndex) = find(SortIndex == BadIndices{FrameIndex}(spotIndex));
end

SortedBestIndices{FrameIndex} = NaN(1,length(BestIndices{FrameIndex}));
for spotIndex = 1:length(BestIndices{FrameIndex})
    SortedBestIndices{FrameIndex}(spotIndex) = find(SortIndex == BestIndices{FrameIndex}(spotIndex));
end

SortedMultiIndices{FrameIndex} = NaN(1,length(MultiIndices{FrameIndex}));
for spotIndex = 1:length(MultiIndices{FrameIndex})
    SortedMultiIndices{FrameIndex}(spotIndex) = find(SortIndex == MultiIndices{FrameIndex}(spotIndex));
end

SortedAggregateIndices{FrameIndex} = NaN(1,length(AggregateIndices{FrameIndex}));
for spotIndex = 1:length(AggregateIndices{FrameIndex})
    SortedAggregateIndices{FrameIndex}(spotIndex) = find(SortIndex == AggregateIndices{FrameIndex}(spotIndex));
end

Spots(FrameIndex).Fits = [Spots(FrameIndex).Fits(SortIndex)];
end

%%

GoodSpots = Spots;
BadSpots= Spots;
MultiSpots = Spots;
BestSpots = Spots;
AggregateSpots = Spots;
for i = 1:length(Spots)
    GoodIndices{i} = unique(SortedGoodIndices{i});
    BadIndices{i} = unique(SortedBadIndices{i});
    MultiIndices{i} = unique(SortedMultiIndices{i});
    BestIndices{i} = unique(SortedBestIndices{i});
    AggregateIndices{i} = unique(SortedAggregateIndices{i});
    GoodSpots(i).Fits = Spots(i).Fits(GoodIndices{i}) ;
    BadSpots(i).Fits = Spots(i).Fits(BadIndices{i}) ;
    MultiSpots(i).Fits = Spots(i).Fits(MultiIndices{i}) ;
    BestSpots(i).Fits = Spots(i).Fits(BestIndices{i}) ;
    AggregateSpots(i).Fits = Spots(i).Fits(AggregateIndices{i}) ;
end


outpath = [liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat'];
save(outpath, 'GoodSpots', 'BadSpots', 'MultiSpots', 'BestSpots', 'AggregateSpots',...
    'GoodIndices', 'BadIndices', 'MultiIndices', 'BestIndices', 'AggregateIndices',...
        'UnsortedGoodIndices','UnsortedBadIndices','UnsortedBestIndices','UnsortedAggregateIndices','UnsortedMultiIndices',...
        'UnsortedSpots');

NumGoodSpots = length(GoodIndices{FrameIndex});

if whos(var2str(Spots)).bytes < 2E9
    save([liveExperiment.resultsFolder,...
        filesep, 'Spots.mat'], 'Spots', '-v6');
else
    save([liveExperiment.resultsFolder,...
        filesep, 'Spots.mat'], 'Spots', '-v7.3', '-nocompression');
end

%%