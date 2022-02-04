function ExtractMembraneFurrowDepths(Prefix, varargin)

% close all
% cleanupObj = onCleanup(@myCleanupFun);

liveExperiment = LiveExperiment(Prefix);

FrameInfo = getFrameInfo(liveExperiment);
load([liveExperiment.resultsFolder, filesep, 'MarkAndFindInfo.Mat'], 'MarkAndFindInfo');

NEmbryos = MarkAndFindInfo.NSeries;
RotatedMembraneFile =[liveExperiment.preFolder, filesep, Prefix, '-Membrane_Rotated.tif'];
RotatedHisFile =[liveExperiment.preFolder, filesep, Prefix, '-His_Rotated.tif'];
CustomRotatedMembraneFile =[liveExperiment.preFolder, filesep, Prefix, '-CustomMembrane_Rotated.tif'];
CustomRotatedHisFile =[liveExperiment.preFolder, filesep, Prefix, '-CustomHis_Rotated.tif'];

MembraneMat = imreadStack2(RotatedMembraneFile, liveExperiment.yDim, liveExperiment.xDim, liveExperiment.nFrames);
HisMat = imreadStack2(RotatedHisFile, liveExperiment.yDim, liveExperiment.xDim, liveExperiment.nFrames);
CustomMembraneMat = imreadStack2(CustomRotatedMembraneFile, liveExperiment.yDim, liveExperiment.xDim, liveExperiment.nFrames);
CustomHisMat = imreadStack2(CustomRotatedHisFile, liveExperiment.yDim, liveExperiment.xDim, liveExperiment.nFrames);
HoldMembraneMat = MembraneMat;
ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;



xSize =size(MembraneMat, 2);
ySize = size(MembraneMat, 1);
PixelSize_um = liveExperiment.pixelSize_um;


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

%%
close all
FurrowSelection = true;
UseHistone = false;
UseCustom = true;
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


MemImage = MembraneMat(yRange,xRange,CurrentEmbryo);


DisplayRange=[min(min(MemImage)),max(max(MemImage))];



Overlay=figure;
set(Overlay,'units', 'normalized', 'position',[0.01, .1, .75, .55]);

overlayAxes = axes(Overlay,'Units', 'normalized', 'Position', [0 0 1 1]);


tb = axtoolbar(overlayAxes);
tb.Visible = 'off';


currentCharacter=1;



% Show the first image
imOverlay = imshow(MemImage,DisplayRange,'Parent',overlayAxes);
hold on
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
    
    set(Overlay,'Name',FigureTitle)
    imOverlay.CData = MemImage;
    
    
    try
        caxis(overlayAxes, DisplayRange);
    end
    if (CompiledEmbryos.Approved(CurrentEmbryo))
        set(Overlay,'Color','g')
    else
        set(Overlay,'Color','r')
    end
    NFurrowPoints=size(CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}, 1);
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
        
        PlotHandle{k} = plot(overlayAxes, CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(k,1)-xlims(1),CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(k,2) -ylims(1),...
            'r+', 'MarkerSize',40);
        LineHandle{k} = plot(overlayAxes,[CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(k,1)-xlims(1),CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(k,3)-xlims(1)],...
            [CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(k,2)-ylims(1),CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(k,4)-ylims(1)],...
            'r','LineStyle','-','Marker','+');
        
        
    end
    
    NCellWidthPoints=size(CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}, 1);
    %refresh ellipses plots by destroying and remaking
    if exist('WidthPlotHandle', 'var')
        cellfun(@delete, WidthPlotHandle);
    end
    if exist('WidthLineHandle', 'var')
        cellfun(@delete, WidthLineHandle);
    end
    
    WidthLineHandle = cell(NCellWidthPoints, 1);
    for k=1:NCellWidthPoints
        
        
        WidthLineHandle{k} = plot(overlayAxes,[CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(k,1)-xlims(1),CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(k,3)-xlims(1)],...
            [CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(k,2)-ylims(1),CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(k,4)-ylims(1)],...
            'g','LineStyle','-','Marker','+');
        
        
    end
    for pl=1:length(APmarkers)
        GuidePlotHandles{pl}.YData =  [0 size(imOverlay.CData,2)];
        GuidePlotHandles{pl}.XData  = [APmarkers(pl), APmarkers(pl)];
    end
    
    
    
    %%
    
    tb = axtoolbar(overlayAxes);
    tb.Visible = 'off';
    
    
    ct=waitforbuttonpress;
    currentCharacter=get(Overlay,'currentcharacter');
    currentMouse=get(overlayAxes,'CurrentPoint');
    
    
    if (ct ~=0)&(currentCharacter=='s')
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
        
        save(CompiledEmbryoPath,'CompiledEmbryos');
        
        
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'normal'))
        currentCharacter=1;
        if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=ySize)&(currentMouse(1,1)<=xSize)
            
            %Add a circle to this location with the median radius of the
            %ellipses found in this frame
            
            %(x, y, a, b, theta, maxcontourvalue, time,
            
            [x,y] = ginput(1);
            if FurrowSelection
                CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}(end+1,:)= ...
                    [currentMouse(1,1)+xlims(1),currentMouse(1,2)+ylims(1),currentMouse(1,1)+xlims(1), y+ylims(1),...
                    0,0,0,0,abs(y-currentMouse(1,2)),abs(y-currentMouse(1,2))*PixelSize_um];
            else
                CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}(end+1,:)= ...
                    [currentMouse(1,1)+xlims(1),(currentMouse(1,2)+y)/2+ylims(1),x+xlims(1), (currentMouse(1,2)+y)/2+ylims(1),...
                    0,0,0,0,abs(x-currentMouse(1,1)),abs(x-currentMouse(1,1))*PixelSize_um];
            end
            
        end
        
        
        
        
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'alt'))
        currentCharacter=1;
        if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=ySize)&(currentMouse(1,1)<=xSize)
            %Find out which ellipses we clicked on so we can delete it
            
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
        end
        
    elseif (ct~=0)&(currentCharacter=='.')&(CurrentEmbryo<find(CompiledEmbryos.Approved,1,'last'))
        
        NewIndex = find(1:NEmbryos > CurrentEmbryo & CompiledEmbryos.Approved, 1);
        if isempty(NewIndex)
            CurrentEmbryo = find(CompiledEmbryos.Approved,1,'last');
        else
            CurrentEmbryo = NewIndex;
        end
        
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
        
        
        MemImage = MembraneMat(yRange,xRange,CurrentEmbryo);
        
        
        DisplayRange=[min(min(MemImage)),max(max(MemImage))];
        
        
        
    elseif (ct~=0)&(currentCharacter==',')&(CurrentEmbryo>find(CompiledEmbryos.Approved,1))
        NewIndex = find(1:NEmbryos < CurrentEmbryo & CompiledEmbryos.Approved, 1,'last');
        if isempty(NewIndex)
            CurrentEmbryo = find(CompiledEmbryos.Approved,1);
        else
            CurrentEmbryo = NewIndex;
        end
        
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
        
        
        MemImage = MembraneMat(yRange,xRange,CurrentEmbryo);
        
        
        DisplayRange=[min(min(MemImage)),max(max(MemImage))];
        
        
    elseif (ct~=0)&(currentCharacter=='j')
        iJump=input('Frame to jump to: ');
        if (floor(iJump)>0)&(iJump<=NEmbryos)
            CurrentEmbryo=iJump;
        else
            disp('Frame out of range.');
        end
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
        if FurrowSelection
            CompiledEmbryos.FurrowMeasurements{CurrentEmbryo}=zeros(0,10);
        else
            CompiledEmbryos.CellWidthMeasurements{CurrentEmbryo}=zeros(0,10);
        end
    elseif (ct~=0)&(currentCharacter=='h')    %Reset the contrast
        UseHistone = ~UseHistone;
        if UseHistone
            if UseCustom
                MembraneMat = CustomHisMat;
            else
                MembraneMat = HisMat;
            end
        else
            if UseCustom
                MembraneMat = CustomMembraneMat;
            else
                MembraneMat = HoldMembraneMat;
            end
        end
        
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
        
        
        MemImage = MembraneMat(yRange,xRange,CurrentEmbryo);
        
        
        DisplayRange=[min(min(MemImage)),max(max(MemImage))];
        
    elseif (ct~=0)&(currentCharacter=='u')    %Reset the contrast
        UseCustom = ~UseCustom;
        if UseHistone
            if UseCustom
                MembraneMat = CustomHisMat;
            else
                MembraneMat = HisMat;
            end
        else
            if UseCustom
                MembraneMat = CustomMembraneMat;
            else
                MembraneMat = HoldMembraneMat;
            end
        end
        
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
        
        
        MemImage = MembraneMat(yRange,xRange,CurrentEmbryo);
        
        
        DisplayRange=[min(min(MemImage)),max(max(MemImage))];
    elseif (ct~=0)&(currentCharacter=='w')    %Reset the contrast
        CompiledEmbryos.Approved(CurrentEmbryo) = false;
    elseif (ct~=0)&(currentCharacter=='q')    %Reset the contrast
        CompiledEmbryos.Approved(CurrentEmbryo) = true;
    end
end
close all
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

save(CompiledEmbryoPath,'CompiledEmbryos');
