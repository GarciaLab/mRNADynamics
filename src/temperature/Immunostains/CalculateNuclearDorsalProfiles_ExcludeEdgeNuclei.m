function CalculateNuclearDorsalProfiles(Prefix,UseCustom, varargin)
%%
ShowPlots = false;
if ~exist('UseCustom', 'var')
    UseCustom = true;
end
liveExperiment= LiveExperiment(Prefix);

CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
load(CompiledEmbryoPath);
FrameInfo = getFrameInfo(liveExperiment);
%%
ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;

Channel1 = liveExperiment.Channel1;
Channel2 = liveExperiment.Channel2;
Channel3 = liveExperiment.Channel3;
Channel4 = liveExperiment.Channel4;
Channel5 = liveExperiment.Channel5;


xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;
PixelSize_um = liveExperiment.pixelSize_um;

%Get the nuclei segmentation data
Ellipses = getEllipses(liveExperiment);
FluoInfoPath = [liveExperiment.resultsFolder, filesep, 'EllipsesFluoInfo.mat'];
load(FluoInfoPath, 'EllipsesFluoInfo');
try
    clrmp = single(hsv(20));
    clrmp = clrmp(randperm(length(clrmp)), :);
catch
    %in case the user doesn't have this colormap, just keep going.
end
Channels = {Channel1, Channel2, Channel3, Channel4, Channel5};

if UseCustom
    RotatedHisFile = [liveExperiment.preFolder, filesep, Prefix, '-CustomHis_Rotated.tif'];
    HisFile = [liveExperiment.preFolder, filesep, Prefix, '-CustomHis.tif'];
else
    RotatedHisFile = [liveExperiment.preFolder, filesep, Prefix, '-His_Rotated.tif'];
    HisFile = [liveExperiment.preFolder, filesep, Prefix, '-His.tif'];
end

rotatedHisMat = imreadStack2(RotatedHisFile, liveExperiment.yDim, liveExperiment.xDim,...
    liveExperiment.nFrames);
hisMat = imreadStack2(HisFile, liveExperiment.yDim, liveExperiment.xDim,...
    liveExperiment.nFrames);



InputChannelIndexes = find(contains(Channels, 'input', 'IgnoreCase', true));
HisChannelIndexes = find(contains(Channels, 'his', 'IgnoreCase', true));
if isempty(InputChannelIndexes)
    warning(['No input channel found. Check correct definition in MovieDatabase.',...
        ' Input channels should use the :input notation.'])
    return;
end
ChannelsToIntegrate = unique([HisChannelIndexes InputChannelIndexes]);

movieMat = getMovieMat(liveExperiment);

numFrames = length(FrameInfo);

%%
xSize = size(hisMat,2);
ySize = size(hisMat,1);
close all
nEmbryos = size(hisMat, 3);
%Get information about the image size
% HisImage=imread([PreProcPath,filesep,Prefix,filesep,D(1).name]);
GoodEmbryos = 1:nEmbryos;

EmptyEllipses = zeros(1, nEmbryos, 'logical');
for CurrentEmbryo = GoodEmbryos
    if isempty(Ellipses{CurrentEmbryo})
        EmptyEllipses(CurrentEmbryo) = true;
    end
end

GoodEmbryos = GoodEmbryos(CompiledEmbryos.Approved & ~EmptyEllipses);
APNarrowbins = 0:0.0125:1;
NarrowProfileNucleiFluoInfo = cell(size(EllipsesFluoInfo, 1), length(Channels));
DorsalNarrowAPProfiles = NaN(nEmbryos, length(APNarrowbins), length(Channels));
APbins = 0:0.025:1;
ProfileNucleiFluoInfo = cell(size(EllipsesFluoInfo, 1), length(Channels));
DorsalAPProfiles = NaN(nEmbryos, length(APbins), length(Channels));
AllDorsalNucleiFluoInfo = cell(size(EllipsesFluoInfo, 1), length(Channels));
AllDorsalNarrowNucleiFluoInfo = cell(size(EllipsesFluoInfo, 1), length(Channels));
DorsalAvgAPProfiles = NaN(nEmbryos, length(APbins), length(Channels));
DorsalAvgNarrowAPProfiles = NaN(nEmbryos, length(APNarrowbins), length(Channels));
DorsalStdAPProfiles = NaN(nEmbryos, length(APbins), length(Channels));
DorsalStdNarrowAPProfiles = NaN(nEmbryos, length(APNarrowbins), length(Channels));
DorsalCountAPProfiles = NaN(nEmbryos, length(APbins), length(Channels));
DorsalCountNarrowAPProfiles = NaN(nEmbryos, length(APNarrowbins), length(Channels));


MiddleMeanDorsalAPProfiles = NaN(nEmbryos, length(APbins), length(Channels));
MiddleMeanDorsalNarrowAPProfiles = NaN(nEmbryos, length(APNarrowbins), length(Channels));
MiddleMeanDorsalAvgAPProfiles = NaN(nEmbryos, length(APbins), length(Channels));
MiddleMeanDorsalAvgNarrowAPProfiles = NaN(nEmbryos, length(APNarrowbins), length(Channels));
MiddleMeanDorsalStdAPProfiles = NaN(nEmbryos, length(APbins), length(Channels));
MiddleMeanDorsalStdNarrowAPProfiles = NaN(nEmbryos, length(APNarrowbins), length(Channels));
MiddleMeanDorsalCountAPProfiles = NaN(nEmbryos, length(APbins), length(Channels));
MiddleMeanDorsalCountNarrowAPProfiles = NaN(nEmbryos, length(APNarrowbins), length(Channels));


%%
for CurrentEmbryo = GoodEmbryos
    
    
    HisImage = hisMat(:,:,CurrentEmbryo);
    HisImage2 = rotatedHisMat(:,:,CurrentEmbryo);
    
    
    
    DisplayRangeHis = [min(min(HisImage)), max(max(HisImage))];
    DisplayRangeHis2 = [min(min(HisImage2)), max(max(HisImage2))];
    if ShowPlots
        close all
        
        FullFigure=figure(1);
        set(FullFigure,'units', 'normalized', 'position',[0.01, .3, .9, .5]);
        
        fullAxes = axes(FullFigure,'Units', 'normalized', 'Position', [0 0 1 1]);
        % Overlay=figure(2);
        % set(Overlay,'units', 'normalized', 'position',[0.01, .1, .45, .65]);
        %
        % overlayAxes = axes(Overlay,'Units', 'normalized', 'Position', [0 0 1 1]);
        %
        %
        % Original=figure(3);
        % set(Original,'units', 'normalized', 'position',[0.51, .1, .45, .65]);
        %
        % originalAxes = axes(Original,'Units', 'normalized', 'Position', [0 0 1 1]);
        
        
        
        
        tb = axtoolbar(fullAxes);
        tb.Visible = 'off';
        % tb2 = axtoolbar(originalAxes);
        % tb2.Visible = 'off';
        % tb3 = axtoolbar(originalAxes);
        % tb3.Visible = 'off';
        
        imFull = imshow(HisImage2,DisplayRangeHis2,'Border','Tight','Parent',fullAxes);
        
        SwitchImageType = false;
        hold(fullAxes,'on')
        
        set(0, 'CurrentFigure', FullFigure)
    end
    %Get the information about the centroids
    [NCentroids,~]=size(Ellipses{CurrentEmbryo});
    
    if ShowPlots
        imFull.CData = HisImage2;
        try
            caxis(fullAxes, DisplayRange);
            
        end
        %refresh ellipses plots by destroying and remaking
        if exist('PlotHandle', 'var')
            cellfun(@delete, PlotHandle);
        end
        
        PlotHandle = cell(NCentroids, 1);
    end
    ellipseFrame = double(Ellipses{CurrentEmbryo});
    EllipsesFluoFrame = double(EllipsesFluoInfo{CurrentEmbryo});
    CoordA = CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,:);
    CoordP = CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,:);
    CoordD = CompiledEmbryos.RotatedCoordDs(CurrentEmbryo,:);
    CoordV = CompiledEmbryos.RotatedCoordVs(CurrentEmbryo,:);
    %bIndex = boundary(ellipseFrame(:,1), ellipseFrame(:,2), 1);
    APLength = CoordP(1)-CoordA(1);
    DVLength = CoordV(2)-CoordD(2);
    DorsalPoints = [];
    
    for k = 1:NCentroids
        if (ellipseFrame(k, 1) >= CoordA(1) + 0.0875*APLength) & ...
                (ellipseFrame(k, 1) <= CoordP(1) - 0.0875*APLength)  & ...
                ellipseFrame(k, 2) <= CoordA(2)
            DorsalPoints(end+1) = k;
        end
    end
    
    dorsalEllipses = ellipseFrame(DorsalPoints,:);
    dorsalEllipses = sortrows(dorsalEllipses, 1);
    
    NFluoCentroids = size(EllipsesFluoFrame, 1);
    DorsalPointsFluo = [];
    for k = 1:NFluoCentroids
        if (EllipsesFluoFrame(k, 1, InputChannelIndexes(1)) >= CoordA(1) + 0.0875*APLength) & ...
                (EllipsesFluoFrame(k, 1, InputChannelIndexes(1)) <= CoordP(1) - 0.0875*APLength)  & ...
                EllipsesFluoFrame(k, 2, InputChannelIndexes(1)) <= CoordA(2)
            DorsalPointsFluo(end+1) = k;
        end
    end
    
    dorsalFluo = EllipsesFluoFrame(DorsalPointsFluo, :,:);
    [out, idx] = sortrows(dorsalFluo(:,:,InputChannelIndexes(1)),1);
    dorsalFluo = dorsalFluo(idx,:,:);
    NuclearDiameter = 2*min(dorsalEllipses(:,3));
    
    EmptyCol = NaN(size(dorsalFluo, 1), 1, size(dorsalFluo, 3));
    dorsalFluo = [dorsalFluo EmptyCol];
    if ShowPlots
        for n=1:size(dorsalFluo, 1)
            colorhash = uint8(mod(round(dorsalFluo(n, 1, InputChannelIndexes(1))+dorsalFluo(n, 2, InputChannelIndexes(1))),20)+1);
            PlotHandle{n} = ellipse(2*dorsalFluo(n, 3,InputChannelIndexes(1)), 2*dorsalFluo(n, 3, InputChannelIndexes(1)),...
                dorsalFluo(n, 4, InputChannelIndexes(1)) * (360/(2*pi)), dorsalFluo(n, 1, InputChannelIndexes(1)),...
                dorsalFluo(n, 2, InputChannelIndexes(1)), clrmp(colorhash,:), 10, fullAxes, 0.5);
        end
    end
    
    MinAPbin = find(round(APbins,5) == 0.1);
    MaxAPbin = find(round(APbins,5) == 0.9);
    DiffBins = diff(APbins);
    
    for chIndex = 2:5
        for i = 1:size(dorsalFluo)
            [minval, minidx] = min(abs(dorsalFluo(i,9:13,chIndex)-dorsalFluo(i, 8, chIndex)));
            if minidx > 1 & minidx < 5
                dorsalFluo(i,14,chIndex) = mean(dorsalFluo(i, (minidx+7):(minidx+9), chIndex));
            end
        end
        
        profileCells = [];
        for binIndex = MinAPbin:MaxAPbin
            cellIndexes = find((dorsalFluo(:,6,chIndex) >= (APbins(binIndex)-DiffBins(1)/2)) & (dorsalFluo(:,6,chIndex) < (APbins(binIndex)+DiffBins(1)/2)));
            if ~isempty(cellIndexes)
                [~, subIndex] = max(dorsalFluo(cellIndexes,8, chIndex));
                singleCellIndex = cellIndexes(subIndex);
                profileCells(end+1) = singleCellIndex;
            end
        end
        
        finalDorsalFluo = dorsalFluo(profileCells,:,:);
        
        nonEdgeIndices = [];
        for i = 1:size(finalDorsalFluo, 1)
            [minval, minidx] = min(abs(finalDorsalFluo(i,9:13,chIndex)-finalDorsalFluo(i, 8, chIndex)));
            if minidx > 1 & minidx < 5
                nonEdgeIndices(end+1) = i;
            end
        end
        finalDorsalFluo = finalDorsalFluo(nonEdgeIndices,:,:);
        ProfileNucleiFluoInfo{CurrentEmbryo, chIndex} = finalDorsalFluo(:,:,chIndex);
        
        for binIndex = MinAPbin:MaxAPbin
            cellIndex = find((finalDorsalFluo(:,6, chIndex)>= (APbins(binIndex)-0.025/2)) & (finalDorsalFluo(:,6, chIndex) < (APbins(binIndex)+0.025/2)));
            if ~isempty(cellIndex)
                DorsalAPProfiles(CurrentEmbryo, binIndex, chIndex) = finalDorsalFluo(cellIndex, 8, chIndex);
            end
        end
        
        
        nonEdgeIndices_dorsalFluo = [];
        for i = 1:size(dorsalFluo, 1)
            [minval, minidx] = min(abs(dorsalFluo(i,9:13,chIndex)-dorsalFluo(i, 8, chIndex)));
            if minidx > 1 & minidx < 5
                nonEdgeIndices_dorsalFluo(end+1) = i;
            end
        end
        
        nonEdgeDorsalFluo = dorsalFluo(nonEdgeIndices_dorsalFluo,:,:);
        AllDorsalNucleiFluoInfo{CurrentEmbryo, chIndex} = nonEdgeDorsalFluo(:,:,chIndex);
        
        for binIndex = MinAPbin:MaxAPbin
            cellIndex = find((nonEdgeDorsalFluo(:,6,chIndex) >= (APbins(binIndex)-DiffBins(1)/2)) & (nonEdgeDorsalFluo(:,6,chIndex) < (APbins(binIndex)+DiffBins(1)/2)));
            if ~isempty(cellIndex)
                DorsalAvgAPProfiles(CurrentEmbryo, binIndex, chIndex) = mean(nonEdgeDorsalFluo(cellIndex, 8,chIndex));
                DorsalStdAPProfiles(CurrentEmbryo, binIndex, chIndex) = std(nonEdgeDorsalFluo(cellIndex, 8,chIndex));
                DorsalCountAPProfiles(CurrentEmbryo, binIndex, chIndex) = length(cellIndex);
            end
        end
        
        MinNarrowAPbin = find(round(APNarrowbins,5) == 0.1);
        MaxNarrowAPbin = find(round(APNarrowbins,5) == 0.9);
        DiffNarrowBins = diff(APNarrowbins);
        % Repeat for NarrowAPbins
        profileCells = [];
        for binIndex = MinNarrowAPbin:MaxNarrowAPbin
            cellIndexes = find((dorsalFluo(:,6,chIndex) >= (APNarrowbins(binIndex)-DiffNarrowBins(1)/2)) & (dorsalFluo(:,6,chIndex) < (APNarrowbins(binIndex)+DiffNarrowBins(1)/2)));
            if ~isempty(cellIndexes)
                [~, subIndex] = max(dorsalFluo(cellIndexes,8,chIndex));
                singleCellIndex = cellIndexes(subIndex);
                profileCells(end+1) = singleCellIndex;
            end
        end
        
        finalDorsalFluo = dorsalFluo(profileCells,:,:);
        nonEdgeIndices = [];
        for i = 1:size(finalDorsalFluo, 1)
            [minval, minidx] = min(abs(finalDorsalFluo(i,9:13,chIndex)-finalDorsalFluo(i, 8, chIndex)));
            if minidx > 1 & minidx < 5
                nonEdgeIndices(end+1) = i;
            end
        end
        finalDorsalFluo = finalDorsalFluo(nonEdgeIndices,:,:);
        ProfileNarrowNucleiFluoInfo{CurrentEmbryo, chIndex} = finalDorsalFluo(:,:,chIndex);
        
        for binIndex =  MinNarrowAPbin:MaxNarrowAPbin
            cellIndex = find((finalDorsalFluo(:,6,chIndex) >= (APNarrowbins(binIndex)-DiffNarrowBins(1)/2)) & (finalDorsalFluo(:,6,chIndex) < (APNarrowbins(binIndex)+DiffNarrowBins(1)/2)));
            if ~isempty(cellIndex)
                DorsalNarrowAPProfiles(CurrentEmbryo, binIndex, chIndex) = finalDorsalFluo(cellIndex, 8,chIndex);
            end
        end
        
        % Repeat for All Cells
        nonEdgeIndices_dorsalFluo = [];
        for i = 1:size(dorsalFluo, 1)
            [minval, minidx] = min(abs(dorsalFluo(i,9:13,chIndex)-dorsalFluo(i, 8, chIndex)));
            if minidx > 1 & minidx < 5
                nonEdgeIndices_dorsalFluo(end+1) = i;
            end
        end
        
        nonEdgeDorsalFluo = dorsalFluo(nonEdgeIndices_dorsalFluo,:,:);
        AllDorsalNarrowNucleiFluoInfo{CurrentEmbryo, chIndex} = nonEdgeDorsalFluo(:,:,chIndex);
        
        for binIndex = MinNarrowAPbin:MaxNarrowAPbin
            cellIndex = find((nonEdgeDorsalFluo(:,6,chIndex) >= (APNarrowbins(binIndex)-DiffNarrowBins(1)/2)) & (nonEdgeDorsalFluo(:,6,chIndex) < (APNarrowbins(binIndex)+DiffNarrowBins(1)/2)));
            if ~isempty(cellIndex) 
                DorsalAvgNarrowAPProfiles(CurrentEmbryo, binIndex, chIndex) = mean(nonEdgeDorsalFluo(cellIndex, 8,chIndex));
                DorsalStdNarrowAPProfiles(CurrentEmbryo, binIndex, chIndex) = std(nonEdgeDorsalFluo(cellIndex, 8,chIndex));
                DorsalCountNarrowAPProfiles(CurrentEmbryo, binIndex, chIndex) = length(cellIndex);
            end
        end
        
        
        
        
        
        %% Add Middle Mean Info
        profileCells = [];
        for binIndex = MinAPbin:MaxAPbin
            cellIndexes = find((dorsalFluo(:,6,chIndex) >= (APbins(binIndex)-DiffBins(1)/2)) & (dorsalFluo(:,6,chIndex) < (APbins(binIndex)+DiffBins(1)/2)));
            if ~isempty(cellIndexes)
                [~, subIndex] = max(dorsalFluo(cellIndexes,14, chIndex));
                singleCellIndex = cellIndexes(subIndex);
                profileCells(end+1) = singleCellIndex;
            end
        end
     
        finalDorsalFluo = dorsalFluo(profileCells,:,:);
        finalDorsalFluo = finalDorsalFluo(~isnan(finalDorsalFluo(:, 14, chIndex)), :,:);
        
        for binIndex = MinAPbin:MaxAPbin
            cellIndex = find((finalDorsalFluo(:,6, chIndex)>= (APbins(binIndex)-0.025/2)) & (finalDorsalFluo(:,6, chIndex) < (APbins(binIndex)+0.025/2)));
            if ~isempty(cellIndex)
                MiddleMeanDorsalAPProfiles(CurrentEmbryo, binIndex, chIndex) = finalDorsalFluo(cellIndex, 14, chIndex);
            end
        end

        nonEdgeDorsalFluo = dorsalFluo(~isnan(dorsalFluo(:,14,chIndex)),:,:);

        for binIndex = MinAPbin:MaxAPbin
            cellIndex = find((nonEdgeDorsalFluo(:,6,chIndex) >= (APbins(binIndex)-DiffBins(1)/2)) & (nonEdgeDorsalFluo(:,6,chIndex) < (APbins(binIndex)+DiffBins(1)/2)));
            if ~isempty(cellIndex)
                MiddleMeanDorsalAvgAPProfiles(CurrentEmbryo, binIndex, chIndex) = mean(nonEdgeDorsalFluo(cellIndex, 14,chIndex));
                MiddleMeanDorsalStdAPProfiles(CurrentEmbryo, binIndex, chIndex) = std(nonEdgeDorsalFluo(cellIndex, 14,chIndex));
                MiddleMeanDorsalCountAPProfiles(CurrentEmbryo, binIndex, chIndex) = length(cellIndex);
            end
        end
        
        MinNarrowAPbin = find(round(APNarrowbins,5) == 0.1);
        MaxNarrowAPbin = find(round(APNarrowbins,5) == 0.9);
        DiffNarrowBins = diff(APNarrowbins);
        % Repeat for NarrowAPbins
        profileCells = [];
        for binIndex = MinNarrowAPbin:MaxNarrowAPbin
            cellIndexes = find((dorsalFluo(:,6,chIndex) >= (APNarrowbins(binIndex)-DiffNarrowBins(1)/2)) & (dorsalFluo(:,6,chIndex) < (APNarrowbins(binIndex)+DiffNarrowBins(1)/2)));
            if ~isempty(cellIndexes)
                [~, subIndex] = max(dorsalFluo(cellIndexes,14,chIndex));
                singleCellIndex = cellIndexes(subIndex);
                profileCells(end+1) = singleCellIndex;
            end
        end
        
        finalDorsalFluo = dorsalFluo(profileCells,:,:);
     
        finalDorsalFluo = finalDorsalFluo(~isnan(finalDorsalFluo(:,14,chIndex)),:,:);
        
        for binIndex =  MinNarrowAPbin:MaxNarrowAPbin
            cellIndex = find((finalDorsalFluo(:,6,chIndex) >= (APNarrowbins(binIndex)-DiffNarrowBins(1)/2)) & (finalDorsalFluo(:,6,chIndex) < (APNarrowbins(binIndex)+DiffNarrowBins(1)/2)));
            if ~isempty(cellIndex)
                MiddleMeanDorsalNarrowAPProfiles(CurrentEmbryo, binIndex, chIndex) = finalDorsalFluo(cellIndex, 14,chIndex);
            end
        end
      
        
        nonEdgeDorsalFluo = dorsalFluo(~isnan(dorsalFluo(:,14,chIndex)),:,:);

        
        for binIndex = MinNarrowAPbin:MaxNarrowAPbin
            cellIndex = find((nonEdgeDorsalFluo(:,6,chIndex) >= (APNarrowbins(binIndex)-DiffNarrowBins(1)/2)) & (nonEdgeDorsalFluo(:,6,chIndex) < (APNarrowbins(binIndex)+DiffNarrowBins(1)/2)));
            if ~isempty(cellIndex)
                MiddleMeanDorsalAvgNarrowAPProfiles(CurrentEmbryo, binIndex, chIndex) = mean(nonEdgeDorsalFluo(cellIndex, 14,chIndex));
                MiddleMeanDorsalStdNarrowAPProfiles(CurrentEmbryo, binIndex, chIndex) = std(nonEdgeDorsalFluo(cellIndex, 14,chIndex));
                MiddleMeanDorsalCountNarrowAPProfiles(CurrentEmbryo, binIndex, chIndex) = length(cellIndex);
            end
        end
        
        
    end
end

%%


DorsalNuclearPath = [liveExperiment.resultsFolder, filesep, 'EdgeCorrectedDorsalNuclearProfiles.mat'];
save(DorsalNuclearPath, 'DorsalAvgAPProfiles', 'DorsalAvgNarrowAPProfiles', 'NarrowProfileNucleiFluoInfo',...
    'DorsalAPProfiles', 'DorsalNarrowAPProfiles','ProfileNucleiFluoInfo',...
    'ProfileNarrowNucleiFluoInfo','AllDorsalNucleiFluoInfo','DorsalStdAPProfiles',...
    'DorsalStdNarrowAPProfiles','DorsalCountAPProfiles','DorsalCountNarrowAPProfiles',...
    'MiddleMeanDorsalAvgAPProfiles', 'MiddleMeanDorsalAvgNarrowAPProfiles', 'MiddleMeanDorsalAPProfiles',...
    'MiddleMeanDorsalNarrowAPProfiles','MiddleMeanDorsalStdAPProfiles','MiddleMeanDorsalStdNarrowAPProfiles',...
    'MiddleMeanDorsalCountAPProfiles','MiddleMeanDorsalCountNarrowAPProfiles');

if ShowPlots
    close all
end
