Prefix='2021-06-10-eVasax5-120mer-GFP-lamin-T22_5C_856x856_LA1_Zoom2_zstep025_45uW';
liveExperiment = LiveExperiment(Prefix);
Spots = getSpots(liveExperiment);

if whos(var2str(Spots)).bytes < 2E9
    save([liveExperiment.resultsFolder, filesep, 'SpotsThresh3000.mat'], 'Spots', '-v6')
else
    save([liveExperiment.resultsFolder,filesep, 'SpotsThresh3000.mat'], 'Spots', '-v7.3', '-nocompression')
end

%%
segment3DNanocages(Prefix, 5000);
SpotsT5000 = fit3DGaussiansToAllSpots(Prefix, '1spot');
%%
if whos(var2str(SpotsT5000)).bytes < 2E9
    save([liveExperiment.resultsFolder, filesep, 'SpotsThresh5000.mat'], 'SpotsT5000', '-v6')
else
    save([liveExperiment.resultsFolder,filesep, 'SpotsThresh5000.mat'], 'SpotsT5000', '-v7.3', '-nocompression')
end
%%
SpotsT3000 = Spots;
if whos(var2str(Spots)).bytes < 2E9
    save([liveExperiment.resultsFolder, filesep, 'Spots.mat'], 'Spots', '-v6')
else
    save([liveExperiment.resultsFolder,filesep, 'Spots.mat'], 'Spots', '-v7.3', '-nocompression')
end

%%
StoredInfo = load([liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat']);
%%
frame_index = 1;
for frame_index = 1:length(SpotsT3000)
SpotsT3000Info(frame_index).IDs =zeros(1, 0, 'uint16');
SpotsT3000Info(frame_index).Rows = zeros(1, 0, 'uint16');
SpotsT3000Info(frame_index).Columns = zeros(1, 0, 'uint16');
SpotsT3000Info(frame_index).Zs = zeros(1, 0, 'uint16');

for sp3_index=1:length(SpotsT3000(frame_index).Fits)
for z_index = 1:SpotsT3000(frame_index).Fits(sp3_index).zCount
NumPixelsForZ = size(SpotsT3000(frame_index).Fits(sp3_index).PixelList{z_index}, 1);
SpotsT3000Info(frame_index).IDs = [SpotsT3000Info(frame_index).IDs sp3_index*ones(1, NumPixelsForZ, 'uint16')];
SpotsT3000Info(frame_index).Zs = [SpotsT3000Info(frame_index).Zs SpotsT3000(frame_index).Fits(sp3_index).z(z_index)*ones(1, NumPixelsForZ, 'uint8')];
SpotsT3000Info(frame_index).Rows = [SpotsT3000Info(frame_index).Rows SpotsT3000(frame_index).Fits(sp3_index).PixelList{z_index}(:,1).'];
SpotsT3000Info(frame_index).Columns = [SpotsT3000Info(frame_index).Columns SpotsT3000(frame_index).Fits(sp3_index).PixelList{z_index}(:,2).'];
end
end
end

%%
for frame_index = 1:length(SpotsT5000)
SpotsT5000Info(frame_index).IDs =zeros(1, 0, 'uint16');
SpotsT5000Info(frame_index).Rows = zeros(1, 0, 'uint16');
SpotsT5000Info(frame_index).Columns = zeros(1, 0, 'uint16');
SpotsT5000Info(frame_index).Zs = zeros(1, 0, 'uint16');
for sp5_index=1:length(SpotsT5000(frame_index).Fits)
for z_index = 1:SpotsT5000(frame_index).Fits(sp5_index).zCount
NumPixelsForZ = size(SpotsT5000(frame_index).Fits(sp5_index).PixelList{z_index}, 1);
SpotsT5000Info(frame_index).IDs = [SpotsT5000Info(frame_index).IDs sp5_index*ones(1, NumPixelsForZ, 'uint16')];
SpotsT5000Info(frame_index).Zs = [SpotsT5000Info(frame_index).Zs SpotsT5000(frame_index).Fits(sp5_index).z(z_index)*ones(1, NumPixelsForZ, 'uint8')];
SpotsT5000Info(frame_index).Rows = [SpotsT5000Info(frame_index).Rows SpotsT5000(frame_index).Fits(sp5_index).PixelList{z_index}(:,1).'];
SpotsT5000Info(frame_index).Columns = [SpotsT5000Info(frame_index).Columns SpotsT5000(frame_index).Fits(sp5_index).PixelList{z_index}(:,2).'];
end
end
end

%%

for frame_index = 1:length(SpotsT3000)
    holdVector = num2cell(zeros(1, length(SpotsT3000(frame_index).Fits), 'uint16'));
    [SpotsT3000(frame_index).Fits(:).SpotT5000Index] = deal(holdVector{:});
    SpotsT3000MatchInfo(frame_index).MatchInfo = zeros(1, length(SpotsT3000(frame_index).Fits), 'uint8');
    SpotsT3000MatchInfo(frame_index).MatchingSpots = cell(1, length(SpotsT3000(frame_index).Fits));
    for sp3_index = 1:length(SpotsT3000(frame_index).Fits)
        spot3AllRows = SpotsT3000Info(frame_index).Rows(SpotsT3000Info(frame_index).IDs == sp3_index);
        spot3AllColumns = SpotsT3000Info(frame_index).Columns(SpotsT3000Info(frame_index).IDs == sp3_index);
        spot3AllZs = SpotsT3000Info(frame_index).Zs(SpotsT3000Info(frame_index).IDs == sp3_index);
        spot3MinZ = min(spot3AllZs);
        spot3MaxZ = max(spot3AllZs);
        spot3MinRow = min(spot3AllRows);
        spot3MaxRow = max(spot3AllRows);
        spot3MinColumn = min(spot3AllColumns);
        spot3MaxColumn = max(spot3AllColumns);
        MatchedConditions = SpotsT5000Info(frame_index).Zs >= spot3MinZ &...
             SpotsT5000Info(frame_index).Zs <= spot3MaxZ &...
             SpotsT5000Info(frame_index).Rows >= spot3MinRow &...
             SpotsT5000Info(frame_index).Rows <= spot3MaxRow &...
             SpotsT5000Info(frame_index).Columns >= spot3MinColumn &...
             SpotsT5000Info(frame_index).Columns <= spot3MaxColumn;
        if sum(MatchedConditions) > 0
            spot5IDSubSet = SpotsT5000Info(frame_index).IDs(MatchedConditions);
            spot5RowsSubSet = SpotsT5000Info(frame_index).Rows(MatchedConditions);
            spot5ColumnsSubSet = SpotsT5000Info(frame_index).Columns(MatchedConditions);
            spot5ZsSubSet = SpotsT5000Info(frame_index).Zs(MatchedConditions);
            SpotsT3000MatchInfo(frame_index).MatchInfo(sp3_index) = length(unique(spot5IDSubSet));
            if SpotsT3000MatchInfo(frame_index).MatchInfo(sp3_index)  == 1
                SpotsT3000MatchInfo(frame_index).MatchingSpots{sp3_index} = unique(spot5IDSubSet);
                SpotsT3000(frame_index).Fits(sp3_index).SpotT5000Index = SpotsT3000MatchInfo(frame_index).MatchingSpots{sp3_index}(1);
            else
                SpotsT3000MatchInfo(frame_index).MatchingSpots{sp3_index} = [];
                zList = unique(spot5ZsSubSet);
                for z_index = 1:length(zList)
                    MatchedConditions2 = MatchedConditions & SpotsT5000Info(frame_index).Zs == zList(z_index);
                    IDSupersubSet = SpotsT5000Info(frame_index).IDs(MatchedConditions2);
                    for id_index = 1:length(IDSupersubSet)
                        if ~ismember(IDSupersubSet(id_index), SpotsT3000MatchInfo(frame_index).MatchingSpots{sp3_index})
                           SpotsT3000MatchInfo(frame_index).MatchingSpots{sp3_index} = [SpotsT3000MatchInfo(frame_index).MatchingSpots{sp3_index} IDSupersubSet(id_index)];
                        end
                    end
                end
                SpotsT3000MatchInfo(frame_index).MatchInfo(sp3_index)  = length( SpotsT3000MatchInfo(frame_index).MatchingSpots{sp3_index} );
            end
        end
        
            
    end
end


%%
GoodIndices = StoredInfo.GoodIndices;
BadIndices = StoredInfo.BadIndices;
MultiIndices = StoredInfo.MultiIndices;
AggregateIndices = StoredInfo.AggregateIndices;
BestIndices = StoredInfo.BestIndices;
GoodSpots = StoredInfo.GoodSpots;
BadSpots = StoredInfo.BadSpots;
MultiSpots = StoredInfo.MultiSpots;
AggregateSpots = StoredInfo.AggregateSpots;
BestSpots = StoredInfo.BestSpots;
NumFrames = length(SpotsT3000);
GoodIndices5000 = cell(1, NumFrames);
BadIndices5000 = cell(1, NumFrames);
MultiIndices5000 = cell(1, NumFrames);
AggregateIndices5000 = cell(1, NumFrames);
BestIndices5000 = cell(1, NumFrames);
for FrameIndex = 1:NumFrames
    for i = 1:length(GoodIndices{FrameIndex})
        if SpotsT3000(FrameIndex).Fits(GoodIndices{FrameIndex}(i)).SpotT5000Index > 0
            GoodIndices5000{FrameIndex} = [GoodIndices5000{FrameIndex}  SpotsT3000(FrameIndex).Fits(GoodIndices{FrameIndex}(i)).SpotT5000Index ];
        end
    end
    for i = 1:length(BadIndices{FrameIndex})
        if SpotsT3000(FrameIndex).Fits(BadIndices{FrameIndex}(i)).SpotT5000Index > 0
            BadIndices5000{FrameIndex} = [BadIndices5000{FrameIndex}  SpotsT3000(FrameIndex).Fits(BadIndices{FrameIndex}(i)).SpotT5000Index ];
        end
    end
    for i = 1:length(BestIndices{FrameIndex})
        if SpotsT3000(FrameIndex).Fits(BestIndices{FrameIndex}(i)).SpotT5000Index > 0
            BestIndices5000{FrameIndex} = [BestIndices5000{FrameIndex}  SpotsT3000(FrameIndex).Fits(BestIndices{FrameIndex}(i)).SpotT5000Index ];
        end
    end
    for i = 1:length(AggregateIndices{FrameIndex})
        if SpotsT3000(FrameIndex).Fits(AggregateIndices{FrameIndex}(i)).SpotT5000Index > 0
            AggregateIndices5000{FrameIndex} = [AggregateIndices5000{FrameIndex}  SpotsT3000(FrameIndex).Fits(AggregateIndices{FrameIndex}(i)).SpotT5000Index ];
        end
    end
    for i = 1:length(MultiIndices{FrameIndex})
        if SpotsT3000(FrameIndex).Fits(MultiIndices{FrameIndex}(i)).SpotT5000Index > 0
            MultiIndices5000{FrameIndex} = [MultiIndices5000{FrameIndex}  SpotsT3000(FrameIndex).Fits(MultiIndices{FrameIndex}(i)).SpotT5000Index ];
        end
    end
end
%%
for i = 1:length(Spots)
    GoodIndices{i} = unique(GoodIndices{i});
    BadIndices{i} = unique(BadIndices{i});
    MultiIndices{i} = unique(MultiIndices{i});
    BestIndices{i} = unique(BestIndices{i});
    AggregateIndices{i} = unique(AggregateIndices{i});
    GoodIndices5000{i} = unique(GoodIndices5000{i});
    BadIndices5000{i} = unique(BadIndices5000{i});
    MultiIndices5000{i} = unique(MultiIndices5000{i});
    BestIndices5000{i} = unique(BestIndices5000{i});
    AggregateIndices5000{i} = unique(AggregateIndices5000{i});
    GoodSpots(i).Fits = Spots(i).Fits(GoodIndices{i}) ;
    BadSpots(i).Fits = Spots(i).Fits(BadIndices{i}) ;
    MultiSpots(i).Fits = Spots(i).Fits(MultiIndices{i}) ;
    BestSpots(i).Fits = Spots(i).Fits(BestIndices{i}) ;
    AggregateSpots(i).Fits = Spots(i).Fits(AggregateIndices{i}) ;
    GoodSpots5000(i).Fits = Spots(i).Fits(GoodIndices5000{i}) ;
    BadSpots5000(i).Fits = Spots(i).Fits(BadIndices5000{i}) ;
    MultiSpots5000(i).Fits = Spots(i).Fits(MultiIndices5000{i}) ;
    BestSpots5000(i).Fits = Spots(i).Fits(BestIndices5000{i}) ;
    AggregateSpots5000(i).Fits = Spots(i).Fits(AggregateIndices5000{i}) ;
end


%%
outpath = [liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat'];
save(outpath, 'GoodSpots', 'BadSpots', 'MultiSpots', 'BestSpots', 'AggregateSpots',...
    'GoodIndices', 'BadIndices', 'MultiIndices', 'BestIndices', 'AggregateIndices',...
    'GoodSpots5000', 'BadSpots5000', 'MultiSpots5000', 'BestSpots5000', 'AggregateSpots5000',...
    'GoodIndices5000', 'BadIndices5000', 'MultiIndices5000', 'BestIndices5000', 'AggregateIndices5000');



