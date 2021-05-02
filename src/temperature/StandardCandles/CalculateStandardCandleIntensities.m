clear all

Prefixes = {'2021-03-24-eVasax5-120mer-GFP-lamin_T25C_Zoom2_zstep025_35uW_01',...
    '2021-03-24-eVasax5-120mer-GFP-lamin_T25C_Zoom2_zstep025_45uW_01',...
    '2021-03-24-eVasax5-120mer-GFP-lamin_T25C_Zoom4_zstep025_35uW_01',...
    '2021-03-24-eVasax5-120mer-GFP-lamin_T25C_Zoom4_zstep025_45uW_01',...
    '2021-03-25-eVasax5-120mer-GFP-lamin_T22_5C_856x856_LA1_Zoom2_zstep025_45uW',...
    '2021-03-25-eVasax5-120mer-GFP-lamin_T22_5C_856x856_LA1_Zoom2_zstep025_60uW',...
    '2021-03-25-eVasax5-120mer-GFP-lamin_T22_5C_856x856_LA1_Zoom4_zstep025_45uW',...
    '2021-03-25-eVasax5-120mer-GFP-lamin_T22_5C_856x856_LA1_Zoom4_zstep025_60uW',...
    '2021-03-26-eVasax5-120mer-GFP-lamin_T20C_856x856_LA1_Zoom2_zstep025_45uW',...
    '2021-03-26-eVasax5-120mer-GFP-lamin_T20C_856x856_LA1_Zoom2_zstep025_60uW',...
    '2021-03-26-eVasax5-120mer-GFP-lamin_T20C_856x856_LA1_Zoom4_zstep025_45uW',...
    '2021-03-26-eVasax5-120mer-GFP-lamin_T20C_856x856_LA1_Zoom4_zstep025_60uW',...
    '2021-03-26-eVasax5-120mer-GFP-lamin_T20C_512x512_LA2_Zoom4_zstep0100_60uW',...
    '2021-03-26-eVasax5-120mer-GFP-lamin_T20C_512x512_LA2_Zoom4_zstep0100_60uW_switch',...
    '2021-03-26-eVasax5-120mer-GFP-lamin_T20C_512x512_LA2_Zoom4_zstep0100_45uW',...
    '2021-03-26-eVasax5-120mer-GFP-lamin_T20C_512x512_LA2_Zoom4_zstep0100_45uW_switch'};
NumSets = length(Prefixes);
xDoG = cell(1, NumSets);
yDoG = cell(1, NumSets);
xFit= cell(1, NumSets);
yFit = cell(1, NumSets);
dogFixedAreaIntensity = cell(1, NumSets);
bwIntensity = cell(1,NumSets);
bwArea = cell(1,NumSets);
GaussianIntensity = cell(1,NumSets);
GaussianInfo_A  = cell(1, NumSets);
GaussianError_A  = cell(1, NumSets);
GaussianInfo_x0  = cell(1, NumSets);
GaussianError_x0  = cell(1, NumSets);
GaussianInfo_y0  = cell(1, NumSets);
GaussianError_y0  = cell(1, NumSets);
GaussianInfo_rho  = cell(1, NumSets);
GaussianError_rho  = cell(1, NumSets);
GaussianInfo_sigx  = cell(1, NumSets);
GaussianError_sigx  = cell(1, NumSets);
GaussianInfo_sigy  = cell(1, NumSets);
GaussianError_sigy  = cell(1, NumSets);
GaussianInfo_offset  = cell(1, NumSets);
GaussianError_offset  = cell(1, NumSets);
GaussianInfo_offx  = cell(1, NumSets);
GaussianError_offx  = cell(1, NumSets);
GaussianInfo_offy  = cell(1, NumSets);
GaussianError_offy  = cell(1, NumSets);
SummedResiduals =  cell(1, NumSets);
z = cell(1,NumSets);
StackInfo = cell(1, NumSets);
intArea = cell(1,NumSets);
discardThis = cell(1,NumSets);
Approved = cell(1,NumSets);
DOGIntensity = cell(1,NumSets);
FixedAreaIntensity = cell(1,NumSets);
bwDiameter = cell(1,NumSets);
bwCircularity = cell(1,NumSets);
bwEccentricity = cell(1,NumSets);
bwMajorAxisLength = cell(1,NumSets);
bwMinorAxisLength = cell(1,NumSets);

nearest_neighbor = cell(1, NumSets);

for j = 1:NumSets
    disp(num2str(j))
    liveExperiment{j} = LiveExperiment(Prefixes{j});
    Spots{j} = getSpots(liveExperiment{j});
    load([liveExperiment{j}.resultsFolder,filesep,'UseSliceInfo.mat'],'UseSliceInfo', 'StackIDs','zIDs');
    StackIDs = StackIDs(UseSliceInfo);
    zIDs = zIDs(UseSliceInfo);
    
    
    
    SummedResiduals{j} = [];
    xDoG{j} =  [];
    yDoG{j} =  [];
    xFit{j} =  [];
    yFit{j} =  [];
    dogFixedAreaIntensity{j} =  [];
    bwIntensity{j} =  [];
    bwArea{j} =  [];
    GaussianIntensity{j} =  [];
    z{j} =  [];
    StackInfo{j} =  [];
    intArea{j} =  [];
    discardThis{j} = [];
    Approved{j} = [];
    DOGIntensity{j} = [];
    FixedAreaIntensity{j}= [];
    
    bwDiameter{j}= [];
    bwCircularity{j}= [];
    bwEccentricity{j}= [];
    bwMajorAxisLength{j}= [];
    bwMinorAxisLength{j}= [];
    
    
    GaussianInfo_A{j}= [];
    GaussianError_A{j}= [];
    GaussianInfo_x0{j}= [];
    GaussianError_x0{j}= [];
    GaussianInfo_y0{j}= [];
    GaussianError_y0{j}= [];
    GaussianInfo_rho{j}= [];
    GaussianError_rho{j}= [];
    GaussianInfo_sigx{j}= [];
    GaussianError_sigx{j}= [];
    GaussianInfo_sigy{j}= [];
    GaussianError_sigy{j}= [];
    GaussianInfo_offset{j}= [];
    GaussianError_offset{j}= [];
    GaussianInfo_offx{j}= [];
    GaussianError_offx{j}= [];
    GaussianInfo_offy{j}= [];
    GaussianError_offy{j}= [];
    
    
    for idx = 1:length(StackIDs)
        i = StackIDs(idx);
        AllZs = [Spots{j}(i).Fits(:).z];
        zMatches = AllZs == zIDs(idx);
        xDoG{j}  = [xDoG{j}  [Spots{j}(i).Fits(zMatches).xDoG]];
        yDoG{j}  = [yDoG{j}  [Spots{j}(i).Fits(zMatches).yDoG]];
        xFit{j}  = [xFit{j}  [Spots{j}(i).Fits(zMatches).xFit]];
        yFit{j}  = [yFit{j}  [Spots{j}(i).Fits(zMatches).yFit]];
        dogFixedAreaIntensity{j}  = [dogFixedAreaIntensity{j}  [Spots{j}(i).Fits(zMatches).dogFixedAreaIntensity]];
        bwIntensity{j}  = [bwIntensity{j}  [Spots{j}(i).Fits(zMatches).bwIntensity]];
        bwArea{j}  = [bwArea{j}  [Spots{j}(i).Fits(zMatches).bwArea]];
        GaussianIntensity{j}  = [GaussianIntensity{j}  [Spots{j}(i).Fits(zMatches).GaussianIntensity]];
        StackInfo{j} = [StackInfo{j} uint8(StackIDs(idx))*ones(1, length([Spots{j}(i).Fits(zMatches).z]), 'uint8')];
        z{j}  = [z{j}  [Spots{j}(i).Fits(zMatches).z]];
        intArea{j}  = [intArea{j}  [Spots{j}(i).Fits(zMatches).intArea]];
        discardThis{j}  = [discardThis{j}  [Spots{j}(i).Fits(zMatches).discardThis]];
        Approved{j}  = [Approved{j}  [Spots{j}(i).Fits(zMatches).Approved]];
        DOGIntensity{j}  = [DOGIntensity{j}  [Spots{j}(i).Fits(zMatches).DOGIntensity]];
        FixedAreaIntensity{j}  = [FixedAreaIntensity{j}  [Spots{j}(i).Fits(zMatches).FixedAreaIntensity]];
        bwDiameter{j}  = [bwDiameter{j}  [Spots{j}(i).Fits(zMatches).bwDiameter]];
        bwCircularity{j}  = [bwCircularity{j}  [Spots{j}(i).Fits(zMatches).bwCircularity]];
        bwEccentricity{j}  = [bwEccentricity{j}  [Spots{j}(i).Fits(zMatches).bwEccentricity]];
        bwMajorAxisLength{j}  = [bwMajorAxisLength{j}  [Spots{j}(i).Fits(zMatches).bwMajorAxisLength]];
        bwMinorAxisLength{j}  = [bwMinorAxisLength{j}  [Spots{j}(i).Fits(zMatches).bwMinorAxisLength]];
        StartIndex = length(GaussianInfo_A{j});
        zIndices = find(AllZs == zIDs(idx));
        NumNewEntries = length(zIndices);
        GaussianInfo_A{j}= [GaussianInfo_A{j}  zeros(1, NumNewEntries, 'single')];
        GaussianError_A{j}= [GaussianError_A{j}  zeros(1, NumNewEntries, 'single')];
        GaussianInfo_x0{j}= [GaussianInfo_x0{j}  zeros(1, NumNewEntries, 'single')];
        GaussianError_x0{j}= [GaussianError_x0{j}  zeros(1, NumNewEntries, 'single')];
        
        GaussianInfo_y0{j}= [GaussianInfo_y0{j}  zeros(1, NumNewEntries, 'single')];
        GaussianError_y0{j}= [GaussianError_y0{j}  zeros(1, NumNewEntries, 'single')];
        GaussianInfo_rho{j}= [GaussianInfo_rho{j}  zeros(1, NumNewEntries, 'single')];
        GaussianError_rho{j}= [GaussianError_rho{j}  zeros(1, NumNewEntries, 'single')];
        GaussianInfo_sigx{j}= [GaussianInfo_sigx{j}  zeros(1, NumNewEntries, 'single')];
        GaussianError_sigx{j}= [GaussianError_sigx{j}  zeros(1, NumNewEntries, 'single')];
        GaussianInfo_sigy{j}= [GaussianInfo_sigy{j}  zeros(1, NumNewEntries, 'single')];
        GaussianError_sigy{j}= [GaussianError_sigy{j}  zeros(1, NumNewEntries, 'single')];
        GaussianInfo_offset{j}= [GaussianInfo_offset{j}  zeros(1, NumNewEntries, 'single')];
        GaussianError_offset{j}= [GaussianError_offset{j}  zeros(1, NumNewEntries, 'single')];
        GaussianInfo_offx{j}= [GaussianInfo_offx{j}  zeros(1, NumNewEntries, 'single')];
        GaussianError_offx{j}= [GaussianError_offx{j}  zeros(1, NumNewEntries, 'single')];
        GaussianInfo_offy{j}= [GaussianInfo_offy{j}  zeros(1, NumNewEntries, 'single')];
        GaussianError_offy{j}= [GaussianError_offy{j}  zeros(1, NumNewEntries, 'single')];
        SummedResiduals{j} = [SummedResiduals{j} zeros(1, NumNewEntries, 'single')];
        for k = 1:NumNewEntries
            GaussianInfo_A{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianInfo.A;
            GaussianError_A{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianError.A;
            GaussianInfo_x0{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianInfo.x0;
            GaussianError_x0{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianError.x0;
            GaussianInfo_y0{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianInfo.y0;
            GaussianError_y0{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianError.y0;
            GaussianInfo_rho{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianInfo.rho;
            GaussianError_rho{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianError.rho;
            GaussianInfo_sigx{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianInfo.sigma_x;
            GaussianError_sigx{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianError.sigma_x;
            GaussianInfo_sigy{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianInfo.sigma_y;
            GaussianError_sigy{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianError.sigma_y;
            GaussianInfo_offset{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianInfo.offset;
            GaussianError_offset{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianError.offset;
            GaussianInfo_offx{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianInfo.offset_x;
            GaussianError_offx{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianError.offset_x;
            GaussianInfo_offy{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianInfo.offset_y;
            GaussianError_offy{j}(StartIndex+k)= Spots{j}(i).Fits(zIndices(k)).GaussianError.offset_y;
            SummedResiduals{j}(StartIndex+k)= sum(sum(abs(Spots{j}(i).Fits(zIndices(k)).GaussianResiduals)));
        end
        
    end
end


for j = 1:NumSets
    nearest_neighbors{j} = zeros(1, length(z{j}), 'single');
    all_stacks = unique(StackInfo{j});
    all_zs = unique(z{j});
    for Stack = 1:length(all_stacks)
        s = all_stacks(Stack);
        for zIndex = 1:length(all_zs)
            zMatch = all_zs(zIndex);
            SpotMatchIndices = find(StackInfo{j} == s & z{j} == zMatch);
            PointLocations = [GaussianInfo_x0{j}(SpotMatchIndices).' GaussianInfo_y0{j}(SpotMatchIndices).'];
            Distances = squareform(pdist(PointLocations));
            Distances(1:length(SpotMatchIndices)+1:end) = max(max(Distances));
            MinDistances = min(Distances);
            nearest_neighbors{j}(SpotMatchIndices) = MinDistances;
        end
    end
end

save('S:/Gabriella/Dropbox\StandardCandles\20210414\SpotInfo.mat')
%%

for j= 1:NumSets
    close all
    figure(1)
    h1 = histogram(bwArea{j});
    [~, AreaIdx] = max(h1.Values);
    AreaMode(j) = (h1.BinEdges(AreaIdx+1)+h1.BinEdges(AreaIdx))/2;
    AreaDelta(j) = std(bwArea{j})/2;
    
    
    figure(2)
    h2 = histogram(bwCircularity{j});
    [~, CircIdx] = max(h2.Values);
    CircMode(j) = (h2.BinEdges(CircIdx+1)+h2.BinEdges(CircIdx))/2;
    CircDelta(j) = .2*CircMode(j);% std(bwCircularity{j});
    
    figure(3)
    h3 = histogram(bwEccentricity{j});
    EccMode = 0.5;
    EccTop = 0.8;
    
    figure(4)
    h4 = histogram(bwMajorAxisLength{j}./bwMinorAxisLength{j});
    [~, AxisIdx] = max(h4.Values);
    AxisMode(j) = (h4.BinEdges(AxisIdx+1)+h4.BinEdges(AxisIdx))/2;
    AxisDelta(j) = 1.5-AxisMode(j);
    
    figure(5)
    h5 = histogram(bwDiameter{j});
    [~, DiamIdx] = max(h5.Values);
    DiamMode(j) = (h5.BinEdges(DiamIdx+1)+h5.BinEdges(DiamIdx))/2;
    DiamDelta(j) = std(bwDiameter{j});
    
    
    figure(6)
    h6 = histogram(bwIntensity{j});
    [~, FluoIdx] = max(h6.Values);
    FluoMode(j) = (h6.BinEdges(FluoIdx+1)+h6.BinEdges(FluoIdx))/2;
    FluoDelta(j) = std(bwIntensity{j});
    
    figure(7)
    h7 = histogram(GaussianIntensity{j});
    [~, FluoIdx] = max(h7.Values);
    GaussianFluoMode(j) = (h7.BinEdges(FluoIdx+1)+h7.BinEdges(FluoIdx))/2;
    GaussianFluoDelta(j) = std(GaussianIntensity{j});
    
    figure(8)
    h8 = histogram(FixedAreaIntensity{j});
    [~, FluoIdx] = max(h8.Values);
    FixedAreaFluoMode(j) = (h8.BinEdges(FluoIdx+1)+h8.BinEdges(FluoIdx))/2;
    FixedAreaFluoDelta(j) = std(FixedAreaIntensity{j});
    
    figure(9)
    h9 = histogram(DOGIntensity{j});
    [~, FluoIdx] = max(h9.Values);
    DOGFluoMode(j) = (h9.BinEdges(FluoIdx+1)+h9.BinEdges(FluoIdx))/2;
    DOGFluoDelta(j) = std(DOGIntensity{j});
    
    figure(10)
    h10 = histogram(dogFixedAreaIntensity{j});
    [~, FluoIdx] = max(h10.Values);
    DOGFixedAreaFluoMode(j) = (h10.BinEdges(FluoIdx+1)+h10.BinEdges(FluoIdx))/2;
    DOGFixedAreaFluoDelta(j) = std(dogFixedAreaIntensity{j});
    
    
    %ValidIndices = (bwArea{j} >= AreaMode(j)-AreaDelta(j)) & (bwArea{j} <= AreaMode(j)+AreaDelta(j))  & ...
    ValidIndices{j} = (bwArea{j} > 1) &  (bwCircularity{j} >= CircMode(j)-CircDelta(j)) & (bwCircularity{j} <= CircMode(j)+CircDelta(j)) & ...
        (bwMajorAxisLength{j}./bwMinorAxisLength{j} >= AxisMode(j)-AxisDelta(j)) &...
        (bwMajorAxisLength{j}./bwMinorAxisLength{j} <= AxisMode(j)+AxisDelta(j)) & ...
        (bwDiameter{j} >= DiamMode(j)-DiamDelta(j)) & (bwDiameter{j} <= DiamMode(j)+DiamDelta(j)) & ...% ...
        (nearest_neighbors{j} >= 1);
    
    
    close all
    figure(1)
    h1 = histogram(bwArea{j});
    hold on
    h1b = histogram(bwArea{j}(ValidIndices{j}));
    
    
    figure(2)
    h2 = histogram(bwCircularity{j});
    hold on
    h2b = histogram(bwCircularity{j}(ValidIndices{j}));
    
    figure(3)
    h3 = histogram(bwEccentricity{j});
    hold on
    h3b = histogram(bwEccentricity{j}(ValidIndices{j}));
    
    figure(4)
    h4 = histogram(bwMajorAxisLength{j}./bwMinorAxisLength{j});
    hold on
    RatioArray = bwMajorAxisLength{j}./bwMinorAxisLength{j};
    h4b = histogram(RatioArray(ValidIndices{j}));
    
    figure(5)
    h5 = histogram(bwDiameter{j});
    hold on
    h5b = histogram(bwDiameter{j}(ValidIndices{j}));
    
    
    figure(6)
    
    h6b = histogram(bwIntensity{j}(ValidIndices{j}));
    [~, FilteredFluoIdx] = max(h6b.Values);
    FilteredFluoMode(j) = (h6b.BinEdges(FilteredFluoIdx+1)+h6b.BinEdges(FilteredFluoIdx))/2;
    FilteredFluoDelta(j) = std(bwIntensity{j}(ValidIndices{j}));
    
    figure(7)
    h7b = histogram(GaussianIntensity{j}(ValidIndices{j}));
    [~, FluoIdx] = max(h7b.Values);
    FilteredGaussianFluoMode(j) = (h7b.BinEdges(FluoIdx+1)+h7b.BinEdges(FluoIdx))/2;
    FiltredGaussianFluoDelta(j) = std(GaussianIntensity{j}(ValidIndices{j}));
    
    figure(8)
    h8b = histogram(FixedAreaIntensity{j}(ValidIndices{j}));
    [~, FluoIdx] = max(h8b.Values);
    FilteredFixedAreaFluoMode(j) = (h8b.BinEdges(FluoIdx+1)+h8b.BinEdges(FluoIdx))/2;
    FilteredFixedAreaFluoDelta(j) = std(FixedAreaIntensity{j}(ValidIndices{j}));
    
    figure(9)
    h9b = histogram(DOGIntensity{j}(ValidIndices{j}));
    [~, FluoIdx] = max(h9b.Values);
    FilteredDOGFluoMode(j) = (h9b.BinEdges(FluoIdx+1)+h9b.BinEdges(FluoIdx))/2;
    FilteredDOGFluoDelta(j) = std(DOGIntensity{j}(ValidIndices{j}));
    
    figure(10)
    h10b = histogram(dogFixedAreaIntensity{j}(ValidIndices{j}));
    [~, FluoIdx] = max(h10b.Values);
    FilteredDOGFixedAreaFluoMode(j) = (h10b.BinEdges(FluoIdx+1)+h10b.BinEdges(FluoIdx))/2;
    FilteredDOGFixedAreaFluoDelta(j) = std(dogFixedAreaIntensity{j}(ValidIndices{j}));
    
    
    disp(Prefixes{j})
    disp('')
end

%%






%%
close all
figure(1)

histogram(dogFixedAreaIntensity{j})

figure(2)
histogram(bwIntensity{j})
%set(gca, 'Yscale', 'Log')


figure(3)
histogram(GaussianIntensity{j})
%set(gca, 'Yscale', 'Log')


figure(4)
histogram(DOGIntensity{j})

figure(5)
histogram(FixedAreaIntensity{j})
%set(gca, 'Yscale', 'Log')




%%
MaxInt = NaN(1,NumSets);
MinInt = NaN(1,NumSets);
MinValue = NaN(1,NumSets);
MaxValue = NaN(1,NumSets);
BinEdges = cell(1,NumSets);

binwidth = 40;
for j=1:NumSets
    MaxInt(j) = ceil(max(DOGIntensity{j}));
    MinInt(j) = floor(min(DOGIntensity{j}));
    MinValue(j) = floor(MinInt(j)/binwidth)*binwidth;
    MaxValue(j) = ceil(MaxInt(j)/binwidth)*binwidth;
    BinEdges{j} = MinValue(j):binwidth:MaxValue(j);
end
close all
figure(1)
h1 = histogram(DOGIntensity{1}, BinEdges{1});
hold on
h2 = histogram(DOGIntensity{2}, BinEdges{2});
hold off




%%
close all
TempAxes = cell(11, 3,4);
plotvalues = cell(11, 3,4);
TempFigure{1} = figure;
set(TempFigure{1} , 'units', 'normalized', 'position', [0.1, 0.1, .8, .8]);
TempAxes{1,1,1} = subplot(3,4,1);
plotvalues{1,1,1} = histogram(GaussianInfo_A{7}-GaussianInfo_offset{7}, 'EdgeColor', 'none');
xlim([0 50])
xlabel('A-offset')
ylabel('Counts')
title('Zoom 4 45 uW 22.5 C')

TempAxes{1,1,2} = subplot(3,4,2);
plotvalues{1,1,2} = histogram(GaussianInfo_A{7}, 'EdgeColor', 'none');
xlim([0 50])
xlabel('A')
ylabel('Counts')
title('Zoom 4 45 uW 22.5 C')

TempAxes{1,1,3} = subplot(3,4,3);
plotvalues{1,1,3} = histogram(GaussianError_A{7}./GaussianInfo_A{7}, 'EdgeColor', 'none', 'BinEdges',...
    0:0.05:max(GaussianError_A{7}./GaussianInfo_A{7}));
xlim([0 1])
xlabel('Fano Factor')
ylabel('Counts')
title('Zoom 4 45 uW 22.5 C')

TempAxes{1,1,4} = subplot(3,4,4);
plotvalues{1,1,4} = scatter(GaussianError_A{7}./GaussianInfo_A{7}, nearest_neighbors{7}, 'r.');
xlim([0 1])
xlabel('Fano Factor')
ylabel('Pixel Distance to Nearest Neighbor')
title('Zoom 4 45 uW 22.5 C')

TempAxes{1,2,1} = subplot(3,4,5);
plotvalues{1,2,1} = histogram(GaussianInfo_A{8}-GaussianInfo_offset{8}, 'EdgeColor', 'none');
xlim([0 50])
xlabel('A-offset')
ylabel('Counts')
title('Zoom 4 60 uW 22.5 C')

TempAxes{1,2,2} = subplot(3,4,6);
plotvalues{1,2,2} = histogram(GaussianInfo_A{8}, 'EdgeColor', 'none');
xlim([0 50])
xlabel('A')
ylabel('Counts')
title('Zoom 4 60 uW 22.5 C')

TempAxes{1,2,3} = subplot(3,4,7);
plotvalues{1,2,3} = histogram(GaussianError_A{8}./GaussianInfo_A{8}, 'EdgeColor', 'none', 'BinEdges',...
    0:0.05:max(GaussianError_A{8}./GaussianInfo_A{8}));
xlim([0 1])
xlabel('Fano Factor')
ylabel('Counts')
title('Zoom 4 60 uW 22.5 C')

TempAxes{1,2,4} = subplot(3,4,8);
plotvalues{1,2,4} = scatter(GaussianError_A{8}./GaussianInfo_A{8}, nearest_neighbors{8}, 'r.');
xlim([0 1])
xlabel('Fano Factor')
ylabel('Pixel Distance to Nearest Neighbor')
title('Zoom 4 60 uW 22.5 C')

TempAxes{1,3,1} = subplot(3,4,9);
plotvalues{1,3,1} = histogram(GaussianInfo_A{7}-GaussianInfo_offset{7}, 'EdgeColor', 'none');
hold on 
histogram(GaussianInfo_A{8}-GaussianInfo_offset{8}, 'EdgeColor', 'none');
hold off
xlim([0 50])
xlabel('A-offset')
ylabel('Counts')


TempAxes{1,3,2} = subplot(3,4,10);
plotvalues{1,3,2} = histogram(GaussianInfo_A{7}, 'EdgeColor', 'none');
hold on 
histogram(GaussianInfo_A{8}, 'EdgeColor', 'none');
hold off
xlim([0 50])
xlabel('A')
ylabel('Counts')


TempAxes{1,3,3} = subplot(3,4,11);
plotvalues{1,3,3} = histogram(GaussianError_A{7}./GaussianInfo_A{7}, 'EdgeColor', 'none', 'BinEdges',...
    0:0.05:max(GaussianError_A{7}./GaussianInfo_A{7}));
hold on 
histogram(GaussianError_A{8}./GaussianInfo_A{8}, 'EdgeColor', 'none', 'BinEdges',...
    0:0.05:max(GaussianError_A{8}./GaussianInfo_A{8}));
hold off
xlim([0 1])
xlabel('Fano Factor')
ylabel('Counts')


TempAxes{1,3,4} = subplot(3,4,12);
plotvalues{1,3,4} = scatter(GaussianInfo_A{7}, nearest_neighbors{7}, 'b.');
hold on 
scatter(GaussianInfo_A{8}, nearest_neighbors{8}, 'r.');
hold off
xlabel('A')
ylabel('Pixel Distance to Nearest Neighbor')

saveas(TempFigure{1} , 'S:/Gabriella/Dropbox\StandardCandles\20210413\Figures\A_22_5C_Plots.png')
%%


TempFigure{2} = figure(2);
set(TempFigure{2} , 'units', 'normalized', 'position', [0.1, 0.1, .8, .8]);
TempAxes{2,1,1} = subplot(2,3,1);
plotvalues{2,1,1} = histogram(SummedResiduals{7}, 'EdgeColor', 'none', 'normalization', 'probability');%, 'BinEdges', 0:.01:max(SummedResiduals{7}));
xlim([0 5000])
xlabel('Summed Residuals')
ylabel('Counts')
title('Zoom 4 45 uW 22.5 C')

TempAxes{2,1,2} = subplot(2,3,2);
plotvalues{2,1,2} = histogram(SummedResiduals{8}, 'EdgeColor', 'none', 'normalization', 'probability');%, 'BinEdges', 0:.01:max(SummedResiduals{8}));
xlim([0 5000])
xlabel('Summed Residuals')
ylabel('Counts')
title('Zoom 4 60 uW 22.5 C')

TempAxes{2,1,3} = subplot(2,3,3);
histogram(SummedResiduals{7}, 'EdgeColor', 'none', 'normalization', 'probability');%, 'BinEdges', 0:.01:max(SummedResiduals{7}));
hold on 
histogram(SummedResiduals{8}, 'EdgeColor', 'none', 'normalization', 'probability');%, 'BinEdges', 0:.01:max(SummedResiduals{8}));
hold off
xlim([0 5000])
xlabel('Summed Residuals')
ylabel('Counts')

TempAxes{2,2,1} = subplot(2,3,4);
scatter(SummedResiduals{7}, GaussianIntensity{7}, 'b.');%, 'BinEdges', 0:.01:max(SummedResiduals{7}));
xlim([0 5000])
xlabel('Summed Residuals')
ylabel('Gaussian Intensity')
title('Zoom 4 45 uW 22.5 C')

TempAxes{2,2,2} = subplot(2,3,5);
scatter(SummedResiduals{8}, GaussianIntensity{8}, 'r.');%, 'BinEdges', 0:.01:max(SummedResiduals{7}));
xlim([0 5000])
xlabel('Summed Residuals')
ylabel('Gaussian Intensity')
title('Zoom 4 60 uW 22.5 C')

TempAxes{2,2,3} = subplot(2,3,6);
scatter(SummedResiduals{7}, GaussianIntensity{7}, 'b.');
hold on 
scatter(SummedResiduals{8}, GaussianIntensity{8}, 'r.');%, 'BinEdges', 0:.01:max(SummedResiduals{7}));
hold off
xlim([0 5000])
xlabel('Summed Residuals')
ylabel('Gaussian Intensity')





TempAxes{1,1,3} = subplot(3,4,3);
plotvalues{1,1,3} = histogram(GaussianError_A{7}./GaussianInfo_A{7}, 'EdgeColor', 'none', 'BinEdges',...
    0:0.05:max(GaussianError_A{7}./GaussianInfo_A{7}));
xlim([0 1])
xlabel('Fano Factor')
ylabel('Counts')
title('Zoom 4 45 uW 22.5 C')

TempAxes{1,1,4} = subplot(3,4,4);
plotvalues{1,1,4} = scatter(GaussianError_A{7}./GaussianInfo_A{7}, nearest_neighbors{7}, 'r.');
xlim([0 1])
xlabel('Fano Factor')
ylabel('Pixel Distance to Nearest Neighbor')
title('Zoom 4 45 uW 22.5 C')

TempAxes{1,2,1} = subplot(3,4,5);
plotvalues{1,2,1} = histogram(GaussianInfo_A{8}-GaussianInfo_offset{8}, 'EdgeColor', 'none');
xlim([0 50])
xlabel('A-offset')
ylabel('Counts')
title('Zoom 4 60 uW 22.5 C')

TempAxes{1,2,2} = subplot(3,4,6);
plotvalues{1,2,2} = histogram(GaussianInfo_A{8}, 'EdgeColor', 'none');
xlim([0 50])
xlabel('A')
ylabel('Counts')
title('Zoom 4 60 uW 22.5 C')

TempAxes{1,2,3} = subplot(3,4,7);
plotvalues{1,2,3} = histogram(GaussianError_A{8}./GaussianInfo_A{8}, 'EdgeColor', 'none', 'BinEdges',...
    0:0.05:max(GaussianError_A{8}./GaussianInfo_A{8}));
xlim([0 1])
xlabel('Fano Factor')
ylabel('Counts')
title('Zoom 4 60 uW 22.5 C')

TempAxes{1,2,4} = subplot(3,4,8);
plotvalues{1,2,4} = scatter(GaussianError_A{8}./GaussianInfo_A{8}, nearest_neighbors{8}, 'r.');
xlim([0 1])
xlabel('Fano Factor')
ylabel('Pixel Distance to Nearest Neighbor')
title('Zoom 4 60 uW 22.5 C')

TempAxes{1,3,1} = subplot(3,4,9);
plotvalues{1,3,1} = histogram(GaussianInfo_A{7}-GaussianInfo_offset{7}, 'EdgeColor', 'none');
hold on 
histogram(GaussianInfo_A{8}-GaussianInfo_offset{8}, 'EdgeColor', 'none');
hold off
xlim([0 50])
xlabel('A-offset')
ylabel('Counts')


TempAxes{1,3,2} = subplot(3,4,10);
plotvalues{1,3,2} = histogram(GaussianInfo_A{7}, 'EdgeColor', 'none');
hold on 
histogram(GaussianInfo_A{8}, 'EdgeColor', 'none');
hold off
xlim([0 50])
xlabel('A')
ylabel('Counts')


TempAxes{1,3,3} = subplot(3,4,11);
plotvalues{1,3,3} = histogram(GaussianError_A{7}./GaussianInfo_A{7}, 'EdgeColor', 'none', 'BinEdges',...
    0:0.05:max(GaussianError_A{7}./GaussianInfo_A{7}));
hold on 
histogram(GaussianError_A{8}./GaussianInfo_A{8}, 'EdgeColor', 'none', 'BinEdges',...
    0:0.05:max(GaussianError_A{8}./GaussianInfo_A{8}));
hold off
xlim([0 1])
xlabel('Fano Factor')
ylabel('Counts')


TempAxes{1,3,4} = subplot(3,4,12);
plotvalues{1,3,4} = scatter(GaussianInfo_A{7}, nearest_neighbors{7}, 'b.');
hold on 
scatter(GaussianInfo_A{8}, nearest_neighbors{8}, 'r.');
hold off
xlabel('A')
ylabel('Pixel Distance to Nearest Neighbor')

%%
close all
TempFigure{3} = figure(3);
set(TempFigure{3} , 'units', 'normalized', 'position', [0.1, 0.1, .8, .8]);
TempAxes{3,1,1} = subplot(3,3,1);
plotvalues{3,1,1} = histogram(GaussianIntensity{7}, 'EdgeColor', 'none', 'BinEdges', 0:10:max(GaussianIntensity{7}),...
    'normalization', 'probability');
xlim([0 800])
xlabel('Gaussian Intensity')
ylabel('Counts')
title('Zoom 4 45 uW 22.5 C')

TempAxes{3,1,2} = subplot(3,3,2);
plotvalues{3,1,2} = scatter(GaussianIntensity{7}, nearest_neighbors{7}, 'b.');
xlim([0 800])
xlabel('Gaussian Intensity')
ylabel('Pixel Distance to Nearest Neighbor')
title('Zoom 4 45 uW 22.5 C')

TempAxes{3,1,3} = subplot(3,3,3);
plotvalues{3,1,3} = scatter(GaussianIntensity{7}, SummedResiduals{7}, 'b.');
xlim([0 800])
xlabel('Gaussian Intensity')
ylabel('Summed Residuals')
title('Zoom 4 45 uW 22.5 C')

TempAxes{1,2,2} = subplot(3,4,6);
plotvalues{1,2,2} = histogram(GaussianInfo_A{8}, 'EdgeColor', 'none');
xlim([0 50])
xlabel('A')
ylabel('Counts')
title('Zoom 4 60 uW 22.5 C')

TempAxes{1,2,3} = subplot(3,4,7);
plotvalues{1,2,3} = histogram(GaussianError_A{8}./GaussianInfo_A{8}, 'EdgeColor', 'none', 'BinEdges',...
    0:0.05:max(GaussianError_A{8}./GaussianInfo_A{8}));
xlim([0 1])
xlabel('Fano Factor')
ylabel('Counts')
title('Zoom 4 60 uW 22.5 C')

TempAxes{1,2,4} = subplot(3,4,8);
plotvalues{1,2,4} = scatter(GaussianError_A{8}./GaussianInfo_A{8}, nearest_neighbors{8}, 'r.');
xlim([0 1])
xlabel('Fano Factor')
ylabel('Pixel Distance to Nearest Neighbor')
title('Zoom 4 60 uW 22.5 C')

TempAxes{1,3,1} = subplot(3,4,9);
plotvalues{1,3,1} = histogram(GaussianInfo_A{7}-GaussianInfo_offset{7}, 'EdgeColor', 'none');
hold on 
histogram(GaussianInfo_A{8}-GaussianInfo_offset{8}, 'EdgeColor', 'none');
hold off
xlim([0 50])
xlabel('A-offset')
ylabel('Counts')


TempAxes{1,3,2} = subplot(3,4,10);
plotvalues{1,3,2} = histogram(GaussianInfo_A{7}, 'EdgeColor', 'none');
hold on 
histogram(GaussianInfo_A{8}, 'EdgeColor', 'none');
hold off
xlim([0 50])
xlabel('A')
ylabel('Counts')


TempAxes{1,3,3} = subplot(3,4,11);
plotvalues{1,3,3} = histogram(GaussianError_A{7}./GaussianInfo_A{7}, 'EdgeColor', 'none', 'BinEdges',...
    0:0.05:max(GaussianError_A{7}./GaussianInfo_A{7}));
hold on 
histogram(GaussianError_A{8}./GaussianInfo_A{8}, 'EdgeColor', 'none', 'BinEdges',...
    0:0.05:max(GaussianError_A{8}./GaussianInfo_A{8}));
hold off
xlim([0 1])
xlabel('Fano Factor')
ylabel('Counts')


TempAxes{1,3,4} = subplot(3,4,12);
plotvalues{1,3,4} = scatter(GaussianInfo_A{7}, nearest_neighbors{7}, 'b.');
hold on 
scatter(GaussianInfo_A{8}, nearest_neighbors{8}, 'r.');
hold off
xlabel('A')
ylabel('Pixel Distance to Nearest Neighbor')

saveas(TempFigure{1} , 'S:/Gabriella/Dropbox\StandardCandles\20210413\Figures\GaussianIntensity_22_5C_Plots.png')


