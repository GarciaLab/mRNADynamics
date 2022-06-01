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
NarrowProfileNucleiFluoInfo = cell(size(EllipsesFluoInfo));
DorsalNarrowAPProfiles = NaN(nEmbryos, length(APNarrowbins)-1, length(Channels));
APbins = 0:0.025:1;
ProfileNucleiFluoInfo = cell(size(EllipsesFluoInfo));
DorsalAPProfiles = NaN(nEmbryos, length(APbins)-1, length(Channels));
AllDorsalNucleiFluoInfo = cell(size(EllipsesFluoInfo));
DorsalAvgAPProfiles = NaN(nEmbryos, length(APbins)-1, length(Channels));
DorsalAvgNarrowAPProfiles = NaN(nEmbryos, length(APNarrowbins)-1, length(Channels));
DorsalStdAPProfiles = NaN(nEmbryos, length(APbins)-1, length(Channels));
DorsalStdNarrowAPProfiles = NaN(nEmbryos, length(APNarrowbins)-1, length(Channels));
DorsalCountAPProfiles = NaN(nEmbryos, length(APbins)-1, length(Channels));
DorsalCountNarrowAPProfiles = NaN(nEmbryos, length(APNarrowbins)-1, length(Channels));
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
        if (ellipseFrame(k, 1) > CoordA(1) + 0.1*APLength) & ...
                (ellipseFrame(k, 1) < CoordP(1) - 0.1*APLength)  & ...
                ellipseFrame(k, 2) < CoordA(2)
            DorsalPoints(end+1) = k;
        end
    end
    
    dorsalEllipses = ellipseFrame(DorsalPoints,:);
    dorsalEllipses = sortrows(dorsalEllipses, 1);
    
    NFluoCentroids = size(EllipsesFluoFrame, 1);
    DorsalPointsFluo = [];
    for k = 1:NFluoCentroids
        if (EllipsesFluoFrame(k, 1, InputChannelIndexes(1)) > CoordA(1) + 0.1*APLength) & ...
                (EllipsesFluoFrame(k, 1, InputChannelIndexes(1)) < CoordP(1) - 0.1*APLength)  & ...
                EllipsesFluoFrame(k, 2, InputChannelIndexes(1)) < CoordA(2)
            DorsalPointsFluo(end+1) = k;
        end
    end
    
    
    
    dorsalFluo = EllipsesFluoFrame(DorsalPointsFluo, :,:);
    [out, idx] = sortrows(dorsalFluo(:,:,InputChannelIndexes(1)),1);
    dorsalFluo = dorsalFluo(idx,:,:);
    NuclearDiameter = 2*min(dorsalEllipses(:,3));
    
    
    
    if ShowPlots
        for k=1:size(dorsalFluo, 1)
            n =k;
            %         PlotHandle{k} = drawellipse('Center',[ellipseFrame(n, 1) ellipseFrame(n, 2)],...
            %             'SemiAxes',[ellipseFrame(n, 3) ellipseFrame(n, 4)], ...
            %             'RotationAngle',ellipseFrame(n, 5) * (360/(2*pi)), 'FaceAlpha', 0,...
            %             'InteractionsAllowed', 'none', 'LabelVisible', 'hover', 'Label', num2str(ellipseFrame(n, 9)));
            colorhash = uint8(mod(round(dorsalFluo(n, 1, InputChannelIndexes(1))+dorsalFluo(n, 2, InputChannelIndexes(1))),20)+1);
            PlotHandle{k} = ellipse(2*dorsalFluo(n, 3,InputChannelIndexes(1)), 2*dorsalFluo(n, 3, InputChannelIndexes(1)),...
                dorsalFluo(n, 4, InputChannelIndexes(1)) * (360/(2*pi)), dorsalFluo(n, 1, InputChannelIndexes(1)),...
                dorsalFluo(n, 2, InputChannelIndexes(1)), clrmp(colorhash,:), 10, fullAxes, 0.5);
            
        end
    end
    
    profileCells = [];
    for binIndex = 1:length(APbins)-1
        cellIndexes = find(dorsalFluo(:,6,InputChannelIndexes(end)) >= APbins(binIndex) & dorsalFluo(:,6,InputChannelIndexes(end)) < APbins(binIndex+1));
        if ~isempty(cellIndexes)
            [~, subIndex] = max(dorsalFluo(cellIndexes,8,InputChannelIndexes(end)));
            singleCellIndex = cellIndexes(subIndex);
            profileCells(end+1) = singleCellIndex;
        end
    end
    
    finalDorsalFluo = dorsalFluo(profileCells,:,:);
    ProfileNucleiFluoInfo{CurrentEmbryo} = finalDorsalFluo;
    
    
    
    for binIndex = 1:length(APbins)-1
        cellIndex = find(finalDorsalFluo(:,6,InputChannelIndexes(end)) >= APbins(binIndex) & finalDorsalFluo(:,6,InputChannelIndexes(end)) < APbins(binIndex+1));
        if ~isempty(cellIndex)
            for channelIndex = ChannelsToIntegrate
                DorsalAPProfiles(CurrentEmbryo, binIndex, channelIndex) = finalDorsalFluo(cellIndex, 8,channelIndex);
            end
        end
    end
    
    
    % Repeat for NarrowAPbins
    profileCells = [];
    for binIndex = 1:length(APNarrowbins)-1
        cellIndexes = find(dorsalFluo(:,6,InputChannelIndexes(end)) >= APNarrowbins(binIndex) & dorsalFluo(:,6,InputChannelIndexes(end)) < APNarrowbins(binIndex+1));
        if ~isempty(cellIndexes)
            [~, subIndex] = max(dorsalFluo(cellIndexes,8,InputChannelIndexes(end)));
            singleCellIndex = cellIndexes(subIndex);
            profileCells(end+1) = singleCellIndex;
        end
    end
    
    finalDorsalFluo = dorsalFluo(profileCells,:,:);
    ProfileNarrowNucleiFluoInfo{CurrentEmbryo} = finalDorsalFluo;
    
    
    
    for binIndex = 1:length(APNarrowbins)-1
        cellIndex = find(finalDorsalFluo(:,6,InputChannelIndexes(end)) >= APNarrowbins(binIndex) & finalDorsalFluo(:,6,InputChannelIndexes(end)) < APNarrowbins(binIndex+1));
        if ~isempty(cellIndex)
            for channelIndex = ChannelsToIntegrate
                DorsalNarrowAPProfiles(CurrentEmbryo, binIndex, channelIndex) = finalDorsalFluo(cellIndex, 8,channelIndex);
            end
        end
    end
    
    % Repeat for All Cells
 
    AllDorsalNucleiFluoInfo{CurrentEmbryo} = dorsalFluo;
    
    for binIndex = 1:length(APbins)-1
        cellIndex = find(dorsalFluo(:,6,InputChannelIndexes(end)) >= APbins(binIndex) & dorsalFluo(:,6,InputChannelIndexes(end)) < APbins(binIndex+1));
        if ~isempty(cellIndex)
            for channelIndex = ChannelsToIntegrate
                DorsalAvgAPProfiles(CurrentEmbryo, binIndex, channelIndex) = mean(dorsalFluo(cellIndex, 8,channelIndex));
                DorsalStdAPProfiles(CurrentEmbryo, binIndex, channelIndex) = std(dorsalFluo(cellIndex, 8,channelIndex));
                DorsalCountAPProfiles(CurrentEmbryo, binIndex, channelIndex) = length(cellIndex);
            end
        end
    end
    
    for binIndex = 1:length(APNarrowbins)-1
        cellIndex = find(dorsalFluo(:,6,InputChannelIndexes(end)) >= APNarrowbins(binIndex) & dorsalFluo(:,6,InputChannelIndexes(end)) < APNarrowbins(binIndex+1));
        if ~isempty(cellIndex)
            for channelIndex = ChannelsToIntegrate
                DorsalAvgNarrowAPProfiles(CurrentEmbryo, binIndex, channelIndex) = mean(dorsalFluo(cellIndex, 8,channelIndex));
                DorsalStdNarrowAPProfiles(CurrentEmbryo, binIndex, channelIndex) = std(dorsalFluo(cellIndex, 8,channelIndex));
                DorsalCountNarrowAPProfiles(CurrentEmbryo, binIndex, channelIndex) = length(cellIndex);
            end
        end
    end
    
    disp('');
    

    
    
end
%%
if ShowPlots
    close all
end

DorsalNuclearPath = [liveExperiment.resultsFolder, filesep, 'DorsalNuclearProfiles.mat'];
save(DorsalNuclearPath, 'DorsalAvgAPProfiles', 'DorsalAvgNarrowAPProfiles', 'NarrowProfileNucleiFluoInfo',...
    'DorsalAPProfiles', 'DorsalNarrowAPProfiles','ProfileNucleiFluoInfo',...
    'ProfileNarrowNucleiFluoInfo','AllDorsalNucleiFluoInfo','DorsalStdAPProfiles',...
    'DorsalStdNarrowAPProfiles','DorsalCountAPProfiles','DorsalCountNarrowAPProfiles');
