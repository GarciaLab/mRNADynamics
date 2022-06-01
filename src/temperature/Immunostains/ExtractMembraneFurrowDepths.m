function ExtractMembraneFurrowDepths(Prefix, UseCustom, varargin)
if ~exist('UseCustom', 'var')
    UseCustom = false;
end

% close all
% cleanupObj = onCleanup(@myCleanupFun);
UseZoomMembrane = false;
liveExperiment = LiveExperiment(Prefix);

FrameInfo = getFrameInfo(liveExperiment);
load([liveExperiment.resultsFolder, filesep, 'MarkAndFindInfo.Mat'], 'MarkAndFindInfo');

NEmbryos = MarkAndFindInfo.NSeries;
MembraneZoomResultsFolder = [liveExperiment.resultsFolder, filesep, 'ZoomMembraneInfo'];


MembraneZoomPixelPath = [MembraneZoomResultsFolder, filesep, 'MembraneZoomPixelSize.mat'];
load(MembraneZoomPixelPath,'PixelSize_um');
MembraneZoomImageSizePath = [MembraneZoomResultsFolder, filesep, 'MembraneZoomImageSize.mat'];
load(MembraneZoomImageSizePath);

MedMembraneFile = [liveExperiment.preFolder, filesep,Prefix, '-MedianZoomMembraneRotated.tif'];
MinMembraneFile = [liveExperiment.preFolder, filesep,Prefix, '-MinZoomMembraneRotated.tif'];
% RotatedMembraneFile =[liveExperiment.preFolder, filesep, Prefix, '-Membrane_Rotated.tif'];
% RotatedHisFile =[liveExperiment.preFolder, filesep, Prefix, '-His_Rotated.tif'];
CustomRotatedMembraneFile =[liveExperiment.preFolder, filesep, Prefix, '-CustomMembrane_Rotated.tif'];
% CustomRotatedHisFile =[liveExperiment.preFolder, filesep, Prefix, '-CustomHis_Rotated.tif'];
% ZoomRotatedMembraneFile =[liveExperiment.preFolder, filesep, Prefix, '-ZoomMembrane_Rotated.tif'];


%HisMat = imreadStack2(RotatedHisFile, liveExperiment.yDim, liveExperiment.xDim, liveExperiment.nFrames);
CustomMembraneMat = imreadStack2(CustomRotatedMembraneFile, liveExperiment.yDim, liveExperiment.xDim, liveExperiment.nFrames);
if UseCustom
    RotatedHisFile = [liveExperiment.preFolder, filesep, Prefix, '-CustomHis_Rotated.tif'];
else
    RotatedHisFile = [liveExperiment.preFolder, filesep, Prefix, '-His_Rotated.tif'];
end
hisMat = imreadStack2(RotatedHisFile, liveExperiment.yDim, liveExperiment.xDim, liveExperiment.nFrames);
HisMat = imresize(hisMat, 2);
HoldMembraneMat = CustomMembraneMat;
if isfile(CustomRotatedMembraneFile) & ~UseCustom
    UseZoomMembrane = true;
    MedZoomMembraneMat = imreadStack2(MedMembraneFile, ySize, xSize, NEmbryos);
    MinZoomMembraneMat = imreadStack2(MinMembraneFile, ySize, xSize, NEmbryos);
    ZoomMembraneMat = MedZoomMembraneMat;
    ZoomPixelSize_um = PixelSize_um;
    HoldMembraneMat = ZoomMembraneMat;
end

%%
ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;



xSize =size(HoldMembraneMat, 2);
ySize = size(HoldMembraneMat, 1);

MemPixelSize_um = PixelSize_um;
if ~UseZoomMembrane
    PixelSize_um = liveExperiment.pixelSize_um;
end


CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
load(CompiledEmbryoPath,'CompiledEmbryos');
Ellipses = getEllipses(liveExperiment);

if ~isfield(CompiledEmbryos, 'nc')
    CompiledEmbryos.nc = NaN(1, NEmbryos);
end

if ~isfield(CompiledEmbryos, 'FurrowMeasurements')
    CompiledEmbryos.FurrowMeasurements = cell(1, NEmbryos);
    for i = 1:NEmbryos
        CompiledEmbryos.FurrowMeasurements{i} = zeros(0,10);
    end
end


if ~isfield(CompiledEmbryos, 'CellWidthMeasurements')
    CompiledEmbryos.CellWidthMeasurements = cell(1, NEmbryos);
    for i = 1:NEmbryos
        CompiledEmbryos.CellWidthMeasurements{i} = zeros(0,10);
    end
end

if ~isfield(CompiledEmbryos, 'ZoomFurrowMeasurements')
    CompiledEmbryos.ZoomFurrowMeasurements = cell(1, NEmbryos);
    for i = 1:NEmbryos
        CompiledEmbryos.ZoomFurrowMeasurements{i} = zeros(0,10);
    end
end


if ~isfield(CompiledEmbryos, 'ZoomCellWidthMeasurements')
    CompiledEmbryos.ZoomCellWidthMeasurements = cell(1, NEmbryos);
    for i = 1:NEmbryos
        CompiledEmbryos.ZoomCellWidthMeasurements{i} = zeros(0,10);
    end
end

%%
close all
FurrowSelection = true;
UseCustom = false;
UseMin = false;
if UseCustom
    MembraneMat = CustomMembraneMat;
    CurrentEmbryo=find(CompiledEmbryos.Approved,1);
    DVLength = abs(CompiledEmbryos.RotatedCoordVs(CurrentEmbryo, 2)-CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2));
    ylims = [round(0.9*CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2)),...
        round(DVLength*1/3+CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2))];
    APLength=abs(CompiledEmbryos.RotatedCoordPs(CurrentEmbryo, 1)-CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1));
    xlims = [round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+APLength*.4),...
        round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+APLength*.6)];
    
    APGrid = 0.4:0.025:0.6;
    APmarkers = APGrid*APLength+CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)-xlims(1);
    
    xRange = xlims(1):xlims(2);
    yRange = ylims(1):ylims(2);
    
end
if UseZoomMembrane
    MembraneMat = ZoomMembraneMat;
    CurrentEmbryo=find(CompiledEmbryos.Approved,1);
    DVLength = abs(CompiledEmbryos.MemRotatedCoordVs(CurrentEmbryo, 2)-CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2));
    ylims = [round(0.9*CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2)),...
        round(DVLength*1/3+CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2))];
    APLength=abs(CompiledEmbryos.MemRotatedCoordPs(CurrentEmbryo, 1)-CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1));
    xlims = [round(CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)+APLength*.3),...
        round(CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)+APLength*.7)];
    
    APGrid = 0.3:0.025:0.7;
    APmarkers = APGrid*APLength+CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)-xlims(1);
    
    xRange = xlims(1):xlims(2);
    yRange = ylims(1):ylims(2);
    
    His_DVLength = abs(CompiledEmbryos.RotatedCoordVs(CurrentEmbryo, 2)-CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2));
    His_ylims = [round(0.9*CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2)),...
        round(His_DVLength*1/3+CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2))];
    His_APLength=abs(CompiledEmbryos.RotatedCoordPs(CurrentEmbryo, 1)-CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1));
    His_xlims = [round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+His_APLength*.3),...
        round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+His_APLength*.7)];
    
    His_APGrid = 0.3:0.025:0.7;
    His_APmarkers = His_APGrid*His_APLength+CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)-His_xlims(1);
    
    His_xRange = His_xlims(1):His_xlims(2);
    His_yRange = His_ylims(1):His_ylims(2);
    
    HisImage = HisMat(yRange, xRange, CurrentEmbryo);
    His_DisplayRange=[min(min(HisImage)),max(max(HisImage))];
    
    
end





MemImage = MembraneMat(yRange,xRange,CurrentEmbryo);
MemImageRGB = cat(3, uint8(MemImage), uint8(MemImage), uint8(MemImage));


DisplayRange=[min(min(MemImage)),max(max(MemImage))];





Overlay=figure;
set(Overlay,'units', 'normalized', 'position',[0.1, .4, .7, .3]);

overlayAxes = axes(Overlay,'Units', 'normalized', 'Position', [0 0 1 1]);


tb = axtoolbar(overlayAxes);
tb.Visible = 'off';

HisFig=figure(2);
set(HisFig,'units', 'normalized', 'position',[0.1, 0.1, .7, .3]);

hisAxes = axes(HisFig,'Units', 'normalized', 'Position', [0 0 1 1]);


tb2 = axtoolbar(hisAxes);
tb2.Visible = 'off';




currentCharacter=1;



% Show the first image
imOverlay = imshow(MemImage,DisplayRange,'Parent',overlayAxes);
imHis = imshow(HisImage,His_DisplayRange,'Parent',hisAxes);
hold(hisAxes, 'on');
HisGuidePlotHandles = cell(1,length(APmarkers));
for pl=1:length(APmarkers)
    HisGuidePlotHandles{pl} = plot(hisAxes, [APmarkers(pl),APmarkers(pl)], [0 size(imHis.CData,2)], 'b-+');
    HisGuidePlotHandles{pl}.Color(4) =0.3;
end


hold(hisAxes, 'off');
hold(overlayAxes, 'on');
% imHis = imshow(HisImageRGB, DisplayRangeHis);
% AlphaMat = squeeze(double((HisImageRGB(:,:,3))));
% level = graythresh((HisImage-1)/256)*256+1;
% AlphaMat = 0.25*double((HisImage > level));
% imHis.AlphaData = AlphaMat;

%%

GuidePlotHandles = cell(1,length(APmarkers));
for pl=1:length(APmarkers)
    GuidePlotHandles{pl} = plot(overlayAxes, [APmarkers(pl),APmarkers(pl)], [0 size(imOverlay.CData,2)], 'b-+');
    GuidePlotHandles{pl}.Color(4) =0.3;
end

%%
try
    if FurrowSelection
        FigureTitle=['Frame: ',num2str(CurrentEmbryo),'/',num2str(NEmbryos),...
            ', nc: ',num2str(CompiledEmbryos.nc(CurrentEmbryo)), ', Add Furrow Depth Measurements'];
    else
        FigureTitle=['Frame: ',num2str(CurrentEmbryo),'/',num2str(NEmbryos),...
            ', nc: ',num2str(CompiledEmbryos.nc(CurrentEmbryo)), ', Add Cell Width Measurements'];
    end
catch
    if FurrowSelection
        FigureTitle=['Frame: ',num2str(CurrentEmbryo),'/',num2str(NEmbryos),...
            ', Add Furrow Depth Measurements'];
    else
        FigureTitle=['Frame: ',num2str(CurrentEmbryo),'/',num2str(NEmbryos),...
            ', Add Cell Width Measurements'];
    end
end
%%
title(overlayAxes, FigureTitle);

hold off


set(0, 'CurrentFigure', Overlay)



%%

while (currentCharacter~='x')
    
    try
        if FurrowSelection
            FigureTitle=['Frame: ',num2str(CurrentEmbryo),'/',num2str(NEmbryos),...
                ', nc: ',num2str(CompiledEmbryos.nc(CurrentEmbryo)), ', Add Furrow Depth Measurements'];
        else
            FigureTitle=['Frame: ',num2str(CurrentEmbryo),'/',num2str(NEmbryos),...
                ', nc: ',num2str(CompiledEmbryos.nc(CurrentEmbryo)), ', Add Cell Width Measurements'];
        end
    catch
        if FurrowSelection
            FigureTitle=['Frame: ',num2str(CurrentEmbryo),'/',num2str(NEmbryos),...
                ', Add Furrow Depth Measurements'];
        else
            FigureTitle=['Frame: ',num2str(CurrentEmbryo),'/',num2str(NEmbryos),...
                ', Add Cell Width Measurements'];
        end
    end
    
    if ~UseZoomMembrane
        NFurrowPoints=size(CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}, 1);
    else
        NFurrowPoints=size(CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}, 1);
    end
    
    hold(hisAxes, 'on')
     %refresh ellipses plots by destroying and remaking
    if exist('HisPlotHandle', 'var')
        cellfun(@delete, HisPlotHandle);
    end
    if exist('HisLineHandle', 'var')
        cellfun(@delete, HisLineHandle);
    end
    HisPlotHandle = cell(NFurrowPoints, 1);
    HisLineHandle = cell(NFurrowPoints, 1);
    for k=1:NFurrowPoints
        if UseZoomMembrane
            HisPlotHandle{k} = plot(hisAxes, CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(k,1),CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(k,2),...
                'r+', 'MarkerSize',40);
            HisLineHandle{k} = plot(hisAxes,[CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(k,1),CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(k,3)],...
                [CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(k,2),CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(k,4)],...
                'r','LineStyle','-','Marker','+');
        end
        
    end
    
    hold(hisAxes, 'off')
    
    set(0, 'CurrentFigure', Overlay)
    hold(overlayAxes, 'on')
    set(Overlay,'Name',FigureTitle)
    imOverlay.CData = MemImage;
    imHis.CData = HisImage;
    
    try
        caxis(overlayAxes, DisplayRange);
    end
    
    try
        caxis(hisAxes, His_DisplayRange);
    end
    
    if (CompiledEmbryos.Approved(CurrentEmbryo))
        set(Overlay,'Color','g')
    else
        set(Overlay,'Color','r')
    end
    
    
    %refresh ellipses plots by destroying and remaking
    if exist('PlotHandle', 'var')
        cellfun(@delete, PlotHandle);
    end
    if exist('LineHandle', 'var')
        cellfun(@delete, LineHandle);
    end
    PlotHandle = cell(NFurrowPoints, 1);
    LineHandle = cell(NFurrowPoints, 1);
    for k=1:NFurrowPoints
        if ~UseZoomMembrane
            PlotHandle{k} = plot(overlayAxes, CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(k,1)-xlims(1),CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(k,2) -ylims(1),...
                'r+', 'MarkerSize',40);
            LineHandle{k} = plot(overlayAxes,[CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(k,1)-xlims(1),CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(k,3)-xlims(1)],...
                [CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(k,2)-ylims(1),CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(k,4)-ylims(1)],...
                'r','LineStyle','-','Marker','+');
        else
            PlotHandle{k} = plot(overlayAxes, CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(k,1),CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(k,2),...
                'r+', 'MarkerSize',40);
            LineHandle{k} = plot(overlayAxes,[CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(k,1),CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(k,3)],...
                [CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(k,2),CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(k,4)],...
                'r','LineStyle','-','Marker','+');
        end
        
    end
    for pl=1:length(APmarkers)
        GuidePlotHandles{pl}.YData =  [0 size(imOverlay.CData,2)];
        GuidePlotHandles{pl}.XData  = [APmarkers(pl), APmarkers(pl)];
    end
    if ~UseZoomMembrane
        
        NCellWidthPoints=size(CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}, 1);
    else
        NCellWidthPoints=size(CompiledEmbryos.ZoomCellWidthMeasurements{CurrentEmbryo}, 1);
    end
    
    
    %refresh ellipses plots by destroying and remaking
    if exist('WidthPlotHandle', 'var')
        cellfun(@delete, WidthPlotHandle);
    end
    if exist('WidthLineHandle', 'var')
        cellfun(@delete, WidthLineHandle);
    end
    
    WidthLineHandle = cell(NCellWidthPoints, 1);
    for k=1:NCellWidthPoints
        
        if ~UseZoomMembrane
            WidthLineHandle{k} = plot(overlayAxes,[CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(k,1)-xlims(1),CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(k,3)-xlims(1)],...
                [CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(k,2)-ylims(1),CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(k,4)-ylims(1)],...
                'g','LineStyle','-','Marker','+');
        else
            WidthLineHandle{k} = plot(overlayAxes,[CompiledEmbryos.ZoomCellWidthMeasurements{CurrentEmbryo}(k,1),CompiledEmbryos.ZoomCellWidthMeasurements{CurrentEmbryo}(k,3)],...
                [CompiledEmbryos.ZoomCellWidthMeasurements{CurrentEmbryo}(k,2),CompiledEmbryos.ZoomCellWidthMeasurements{CurrentEmbryo}(k,4)],...
                'g','LineStyle','-','Marker','+');
        end
        
        
    end
    hold(overlayAxes, 'off')
    
    %%
    
    tb = axtoolbar(overlayAxes);
    tb.Visible = 'off';
    
    
    ct=waitforbuttonpress;
    currentCharacter=get(Overlay,'currentcharacter');
    currentMouse=get(overlayAxes,'CurrentPoint');
    
    
    if (ct ~=0)&(currentCharacter=='s')
        if ~UseZoomMembrane
            for embryo_index = 1:NEmbryos
                if ~isempty(CompiledEmbryos.FurrowMeasurements{CurrentEmbryo})
                    CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(:,5:6) = ...
                        TransformToOriginalImageCoordinates(CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(:,1:2), Prefix, CompiledEmbryos, liveExperiment, embryo_index);
                    CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(:,7:8) = ...
                        TransformToOriginalImageCoordinates(CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(:,3:4), Prefix, CompiledEmbryos, liveExperiment, embryo_index);
                end
                if ~isempty(CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo})
                    CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(:,5:6) = ...
                        TransformToOriginalImageCoordinates(CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(:,1:2), Prefix, CompiledEmbryos, liveExperiment, embryo_index);
                    CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(:,7:8) = ...
                        TransformToOriginalImageCoordinates(CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(:,3:4), Prefix, CompiledEmbryos, liveExperiment, embryo_index);
                end
            end
        end
        save(CompiledEmbryoPath,'CompiledEmbryos');
        
        
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'normal'))
        currentCharacter=1;
        if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=ySize)&(currentMouse(1,1)<=xSize)
            
            %Add a circle to this location with the median radius of the
            %ellipses found in this frame
            
            %(x, y, a, b, theta, maxcontourvalue, time,
            
            [x,y] = ginput(1);
            if ~UseZoomMembrane
                if FurrowSelection
                    CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(end+1,:)= ...
                        [currentMouse(1,1)+xlims(1),currentMouse(1,2)+ylims(1),currentMouse(1,1)+xlims(1), y+ylims(1),...
                        0,0,0,0,abs(y-currentMouse(1,2)),abs(y-currentMouse(1,2))*PixelSize_um];
                else
                    CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(end+1,:)= ...
                        [currentMouse(1,1)+xlims(1),(currentMouse(1,2)+y)/2+ylims(1),x+xlims(1), (currentMouse(1,2)+y)/2+ylims(1),...
                        0,0,0,0,abs(x-currentMouse(1,1)),abs(x-currentMouse(1,1))*PixelSize_um];
                end
            else
                if FurrowSelection
                    CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(end+1,:)= ...
                        [currentMouse(1,1),currentMouse(1,2),currentMouse(1,1), y,...
                        0,0,0,0,abs(y-currentMouse(1,2)),abs(y-currentMouse(1,2))*MemPixelSize_um];
                else
                    CompiledEmbryos.ZoomCellWidthMeasurements{CurrentEmbryo}(end+1,:)= ...
                        [currentMouse(1,1),(currentMouse(1,2)+y)/2,x, (currentMouse(1,2)+y)/2,...
                        0,0,0,0,abs(x-currentMouse(1,1)),abs(x-currentMouse(1,1))*MemPixelSize_um];
                end
            end
            
        end
        
        
        
        
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'alt'))
        currentCharacter=1;
        if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=ySize)&(currentMouse(1,1)<=xSize)
            %Find out which ellipses we clicked on so we can delete it
            
            if ~UseZoomMembrane
                if FurrowSelection
                    %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
                    Distances=sqrt((CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(:,1)-(currentMouse(1,1)+xlims(1))).^2+...
                        (CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(:,2)-(currentMouse(1,2)+ylims(1))).^2);
                    if ~isempty(Distances)
                        [~,MinIndex]=min(Distances);
                        
                        CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}=[CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(1:MinIndex-1,:);...
                            CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(MinIndex+1:end,:)];
                    end
                else
                    Distances=sqrt((CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(:,1)-(currentMouse(1,1)+xlims(1))).^2+...
                        (CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(:,2)-(currentMouse(1,2)+ylims(1))).^2);
                    if ~isempty(Distances)
                        [~,MinIndex]=min(Distances);
                        
                        CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}=[CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(1:MinIndex-1,:);...
                            CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(MinIndex+1:end,:)];
                    end
                    
                end
            else
                if FurrowSelection
                    %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
                    Distances=sqrt((CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(:,1)-(currentMouse(1,1))).^2+...
                        (CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(:,2)-(currentMouse(1,2))).^2);
                    if ~isempty(Distances)
                        [~,MinIndex]=min(Distances);
                        
                        CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}=[CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(1:MinIndex-1,:);...
                            CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}(MinIndex+1:end,:)];
                    end
                else
                    Distances=sqrt((CompiledEmbryos.ZoomCellWidthMeasurements{CurrentEmbryo}(:,1)-(currentMouse(1,1))).^2+...
                        (CompiledEmbryos.ZoomCellWidthMeasurements{CurrentEmbryo}(:,2)-(currentMouse(1,2))).^2);
                    if ~isempty(Distances)
                        [~,MinIndex]=min(Distances);
                        
                        CompiledEmbryos.ZoomCellWidthMeasurements{CurrentEmbryo}=[CompiledEmbryos.ZoomCellWidthMeasurements{CurrentEmbryo}(1:MinIndex-1,:);...
                            CompiledEmbryos.ZoomCellWidthMeasurements{CurrentEmbryo}(MinIndex+1:end,:)];
                    end
                    
                end
            end
        end
        
    elseif (ct~=0)&(currentCharacter=='.')&(CurrentEmbryo<find(CompiledEmbryos.Approved,1,'last'))
        
        NewIndex = find(1:NEmbryos > CurrentEmbryo & CompiledEmbryos.Approved, 1);
        if isempty(NewIndex)
            CurrentEmbryo = find(CompiledEmbryos.Approved,1,'last');
        else
            CurrentEmbryo = NewIndex;
        end
        
        if UseCustom
            DVLength = abs(CompiledEmbryos.RotatedCoordVs(CurrentEmbryo, 2)-CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2));
            ylims = [round(0.9*CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2)),...
                round(DVLength*1/3+CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2))];
            APLength=abs(CompiledEmbryos.RotatedCoordPs(CurrentEmbryo, 1)-CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1));
            xlims = [round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+APLength*.4),...
                round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+APLength*.6)];
            
            APGrid = 0.4:0.025:0.6;
            APmarkers = APGrid*APLength+CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)-xlims(1);
            
            xRange = xlims(1):xlims(2);
            yRange = ylims(1):ylims(2);
            
        end
        if UseZoomMembrane
            DVLength = abs(CompiledEmbryos.MemRotatedCoordVs(CurrentEmbryo, 2)-CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2));
            ylims = [round(0.9*CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2)),...
                round(DVLength*1/3+CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2))];
            APLength=abs(CompiledEmbryos.MemRotatedCoordPs(CurrentEmbryo, 1)-CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1));
            xlims = [round(CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)+APLength*.3),...
                round(CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)+APLength*.7)];
            
            APGrid = 0.3:0.025:0.7;
            APmarkers = APGrid*APLength+CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)-xlims(1);
            
            xRange = xlims(1):xlims(2);
            yRange = ylims(1):ylims(2);
            
            His_DVLength = abs(CompiledEmbryos.RotatedCoordVs(CurrentEmbryo, 2)-CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2));
    His_ylims = [round(0.9*CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2)),...
        round(His_DVLength*1/3+CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2))];
    His_APLength=abs(CompiledEmbryos.RotatedCoordPs(CurrentEmbryo, 1)-CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1));
    His_xlims = [round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+His_APLength*.3),...
        round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+His_APLength*.7)];
    
    His_APGrid = 0.3:0.025:0.7;
    His_APmarkers = His_APGrid*His_APLength+CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)-His_xlims(1);
    
    His_xRange = His_xlims(1):His_xlims(2);
    His_yRange = His_ylims(1):His_ylims(2);
    
    HisImage = HisMat(yRange, xRange, CurrentEmbryo);
    His_DisplayRange=[min(min(HisImage)),max(max(HisImage))];
            
        end
        
        
        
        
        
        
        MemImage = MembraneMat(yRange,xRange,CurrentEmbryo);
        DisplayRange=[min(min(MemImage)),max(max(MemImage))];
        
    elseif (ct~=0)&(currentCharacter==',')&(CurrentEmbryo>find(CompiledEmbryos.Approved,1))
        NewIndex = find(1:NEmbryos < CurrentEmbryo & CompiledEmbryos.Approved, 1,'last');
        if isempty(NewIndex)
            CurrentEmbryo = find(CompiledEmbryos.Approved,1);
        else
            CurrentEmbryo = NewIndex;
        end
        
        if UseCustom
            
            DVLength = abs(CompiledEmbryos.RotatedCoordVs(CurrentEmbryo, 2)-CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2));
            ylims = [round(0.9*CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2)),...
                round(DVLength*1/3+CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2))];
            APLength=abs(CompiledEmbryos.RotatedCoordPs(CurrentEmbryo, 1)-CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1));
            xlims = [round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+APLength*.4),...
                round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+APLength*.6)];
            
            APGrid = 0.4:0.025:0.6;
            APmarkers = APGrid*APLength+CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)-xlims(1);
            
            xRange = xlims(1):xlims(2);
            yRange = ylims(1):ylims(2);
            
        end
        if UseZoomMembrane
            DVLength = abs(CompiledEmbryos.MemRotatedCoordVs(CurrentEmbryo, 2)-CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2));
            ylims = [round(0.9*CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2)),...
                round(DVLength*1/3+CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2))];
            APLength=abs(CompiledEmbryos.MemRotatedCoordPs(CurrentEmbryo, 1)-CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1));
            xlims = [round(CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)+APLength*.3),...
                round(CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)+APLength*.7)];
            
            APGrid = 0.3:0.025:0.7;
            APmarkers = APGrid*APLength+CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)-xlims(1);
            
            xRange = xlims(1):xlims(2);
            yRange = ylims(1):ylims(2);
            
            His_DVLength = abs(CompiledEmbryos.RotatedCoordVs(CurrentEmbryo, 2)-CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2));
    His_ylims = [round(0.9*CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2)),...
        round(His_DVLength*1/3+CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2))];
    His_APLength=abs(CompiledEmbryos.RotatedCoordPs(CurrentEmbryo, 1)-CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1));
    His_xlims = [round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+His_APLength*.3),...
        round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+His_APLength*.7)];
    
    His_APGrid = 0.3:0.025:0.7;
    His_APmarkers = His_APGrid*His_APLength+CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)-His_xlims(1);
    
    His_xRange = His_xlims(1):His_xlims(2);
    His_yRange = His_ylims(1):His_ylims(2);
    
    HisImage = HisMat(yRange, xRange, CurrentEmbryo);
    His_DisplayRange=[min(min(HisImage)),max(max(HisImage))];
            
            
        end
        
        
        
        
        
        MemImage = MembraneMat(yRange,xRange,CurrentEmbryo);
        DisplayRange=[min(min(MemImage)),max(max(MemImage))];
        
        
        
    elseif (ct~=0)&(currentCharacter=='j')
        iJump=input('Frame to jump to: ');
        if (floor(iJump)>0)&(iJump<=NEmbryos)
            CurrentEmbryo=iJump;
        else
            disp('Frame out of range.');
        end
        if UseCustom
            DVLength = abs(CompiledEmbryos.RotatedCoordVs(CurrentEmbryo, 2)-CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2));
            ylims = [round(0.9*CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2)),...
                round(DVLength*1/3+CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2))];
            APLength=abs(CompiledEmbryos.RotatedCoordPs(CurrentEmbryo, 1)-CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1));
            xlims = [round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+APLength*.4),...
                round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+APLength*.6)];
            
            APGrid = 0.4:0.025:0.6;
            APmarkers = APGrid*APLength+CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)-xlims(1);
            
            xRange = xlims(1):xlims(2);
            yRange = ylims(1):ylims(2);
            
        end
        if UseZoomMembrane
            DVLength = abs(CompiledEmbryos.MemRotatedCoordVs(CurrentEmbryo, 2)-CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2));
            ylims = [round(0.9*CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2)),...
                round(DVLength*1/3+CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2))];
            APLength=abs(CompiledEmbryos.MemRotatedCoordPs(CurrentEmbryo, 1)-CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1));
            xlims = [round(CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)+APLength*.3),...
                round(CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)+APLength*.7)];
            
            APGrid = 0.3:0.025:0.7;
            APmarkers = APGrid*APLength+CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)-xlims(1);
            
            xRange = xlims(1):xlims(2);
            yRange = ylims(1):ylims(2);
            
            His_DVLength = abs(CompiledEmbryos.RotatedCoordVs(CurrentEmbryo, 2)-CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2));
    His_ylims = [round(0.9*CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2)),...
        round(His_DVLength*1/3+CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2))];
    His_APLength=abs(CompiledEmbryos.RotatedCoordPs(CurrentEmbryo, 1)-CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1));
    His_xlims = [round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+His_APLength*.3),...
        round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+His_APLength*.7)];
    
    His_APGrid = 0.3:0.025:0.7;
    His_APmarkers = His_APGrid*His_APLength+CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)-His_xlims(1);
    
    His_xRange = His_xlims(1):His_xlims(2);
    His_yRange = His_ylims(1):His_ylims(2);
    
    HisImage = HisMat(yRange, xRange, CurrentEmbryo);
    His_DisplayRange=[min(min(HisImage)),max(max(HisImage))];
            
            
        end
        
        
        
        
        
        MemImage = MembraneMat(yRange,xRange,CurrentEmbryo);
        DisplayRange=[min(min(MemImage)),max(max(MemImage))];
        
        
    elseif (ct~=0)&(currentCharacter=='4')
        CompiledEmbryos.nc(CurrentEmbryo) = 14;
    elseif (ct~=0)&(currentCharacter=='3')
        CompiledEmbryos.nc(CurrentEmbryo) = 13;
    elseif (ct~=0)&(currentCharacter=='0')
        CompiledEmbryos.nc(CurrentEmbryo) = NaN;
    elseif (ct~=0)&(currentCharacter=='m')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)/1.5;
        
    elseif (ct~=0)&(currentCharacter=='n')    %Decrease contrast
        DisplayRange(2)=DisplayRange(2)*1.5;
        
    elseif (ct~=0)&(currentCharacter=='r')    %Reset the contrast
        DisplayRange=[min(min(MemImage)),max(max(MemImage))];
    elseif (ct~=0)&(currentCharacter=='f')    %Reset the contrast
        FurrowSelection = ~FurrowSelection;
        
    elseif (ct~=0)&(currentCharacter=='d')    %Delete all ellipses in the current frame
        if ~UseZoomMembrane
            if FurrowSelection
                CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}=zeros(0,10);
            else
                CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}=zeros(0,10);
            end
        else
            if FurrowSelection
                CompiledEmbryos.ZoomFurrowMeasurements{CurrentEmbryo}=zeros(0,10);
            else
                CompiledEmbryos.ZoomCellWidthMeasurements{CurrentEmbryo}=zeros(0,10);
            end
        end
        
        
        
    elseif (ct~=0)&(currentCharacter=='w')    %Reset the contrast
        CompiledEmbryos.Approved(CurrentEmbryo) = false;
    elseif (ct~=0)&(currentCharacter=='q')    %Reset the contrast
        CompiledEmbryos.Approved(CurrentEmbryo) = true;
    elseif (ct~=0)&(currentCharacter=='t')    %Switch mat to use
        if UseZoomMembrane
            
            if UseMin
                UseMin = false;
                MembraneMat = MinZoomMembraneMat;
            else
                UseMin = true;
                MembraneMat = MedZoomMembraneMat;
            end
            DVLength = abs(CompiledEmbryos.MemRotatedCoordVs(CurrentEmbryo, 2)-CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2));
            ylims = [round(0.9*CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2)),...
                round(DVLength*1/3+CompiledEmbryos.MemRotatedCoordDs(CurrentEmbryo, 2))];
            APLength=abs(CompiledEmbryos.MemRotatedCoordPs(CurrentEmbryo, 1)-CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1));
            xlims = [round(CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)+APLength*.3),...
                round(CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)+APLength*.7)];
            
            APGrid = 0.3:0.025:0.7;
            APmarkers = APGrid*APLength+CompiledEmbryos.MemRotatedCoordAs(CurrentEmbryo, 1)-xlims(1);
            
            xRange = xlims(1):xlims(2);
            yRange = ylims(1):ylims(2);
            
            
            MemImage = MembraneMat(yRange,xRange,CurrentEmbryo);
            DisplayRange=[min(min(MemImage)),max(max(MemImage))];
            
            His_DVLength = abs(CompiledEmbryos.RotatedCoordVs(CurrentEmbryo, 2)-CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2));
    His_ylims = [round(0.9*CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2)),...
        round(His_DVLength*1/3+CompiledEmbryos.RotatedCoordDs(CurrentEmbryo, 2))];
    His_APLength=abs(CompiledEmbryos.RotatedCoordPs(CurrentEmbryo, 1)-CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1));
    His_xlims = [round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+His_APLength*.3),...
        round(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)+His_APLength*.7)];
    
    His_APGrid = 0.3:0.025:0.7;
    His_APmarkers = His_APGrid*His_APLength+CompiledEmbryos.RotatedCoordAs(CurrentEmbryo, 1)-His_xlims(1);
    
    His_xRange = His_xlims(1):His_xlims(2);
    His_yRange = His_ylims(1):His_ylims(2);
    
    HisImage = HisMat(yRange, xRange, CurrentEmbryo);
    His_DisplayRange=[min(min(HisImage)),max(max(HisImage))];
            
        end
        
        
        
        
        
    end
end
close all
if ~UseZoomMembrane
    for embryo_index = 1:NEmbryos
        if ~isempty(CompiledEmbryos.FurrowMeasurements{CurrentEmbryo})
            CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(:,5:6) = ...
                TransformToOriginalImageCoordinates(CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(:,1:2), Prefix, CompiledEmbryos, liveExperiment, embryo_index);
            CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(:,7:8) = ...
                TransformToOriginalImageCoordinates(CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(:,3:4), Prefix, CompiledEmbryos, liveExperiment, embryo_index);
        end
        if ~isempty(CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo})
            CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(:,5:6) = ...
                TransformToOriginalImageCoordinates(CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(:,1:2), Prefix, CompiledEmbryos, liveExperiment, embryo_index);
            CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(:,7:8) = ...
                TransformToOriginalImageCoordinates(CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(:,3:4), Prefix, CompiledEmbryos, liveExperiment, embryo_index);
        end
    end
end

save(CompiledEmbryoPath,'CompiledEmbryos');
