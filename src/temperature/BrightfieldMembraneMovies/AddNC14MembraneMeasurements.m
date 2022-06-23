function AddNC14MembraneMeasurements(Prefix, reset_values)
close all
if ~exist('reset_values','var')
    
    reset_values =false;
end


liveExperiment = LiveExperiment(Prefix);

FrameInfo = getFrameInfo(liveExperiment);
FrameTimes_seconds = [FrameInfo(:).Time];
FrameTimes_minutes = FrameTimes_seconds/60;

ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;

Channel1 = liveExperiment.Channel1;
anaphaseFrames = liveExperiment.anaphaseFrames;
nc9 = anaphaseFrames(1);
nc10 = anaphaseFrames(2);
nc11 = anaphaseFrames(3);
nc12 = anaphaseFrames(4);
nc13 = anaphaseFrames(5);
nc14 = anaphaseFrames(6);

FrameTimes_NC14 = FrameTimes_minutes(nc14:end)-FrameTimes_minutes(nc14);

xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;
PixelSize_um = liveExperiment.pixelSize_um;

load([DropboxFolder,filesep,Prefix,filesep,'EmbryoOrientationInfo.mat'])


membraneMat = imreadStack2([liveExperiment.preFolder, filesep,...
    liveExperiment.Prefix, '-Membrane_Rotated.tif'], liveExperiment.yDim, liveExperiment.xDim, liveExperiment.nFrames);
membraneMat = double(membraneMat);

nFrames = size(membraneMat, 3);
%Get information about the image size
% HisImage=imread([PreProcPath,filesep,Prefix,filesep,D(1).name]);

if ~reset_values
    if isfile([DropboxFolder,filesep,Prefix,filesep,'FurrowCanalDepthMeasurements.mat'])
        load([DropboxFolder,filesep,Prefix,filesep,'FurrowCanalDepthMeasurements.mat'])
    end
end

if ~exist('FurrowMeasurements','var')
    
    FurrowMeasurements = cell(1,nFrames);
end
if ~exist('PerivitellineMeasurements','var')
    
    PerivitellineMeasurements = cell(1,nFrames);
end
if ~exist('NumFurrowMeasurements','var')
    NumFurrowMeasuremnts= zeros(1,nFrames);
    for frame_index=1:nFrames
        FurrowMeasurements{frame_index} =NaN(0,3,'double');
    end
end
if ~exist('NumPVMeasurements','var')
    NumPVMeasurements= zeros(1,nFrames);

end

if ~exist('DeltaFC_pixels', 'var')
    DeltaFC_pixels = NaN(nFrames,2);
end

if ~exist('DeltaFC_um', 'var')
    DeltaFC_um = NaN(nFrames,2);
end

if ~exist('EmbryoBoundaryPositions', 'var')
    EmbryoBoundaryPositions = BoundaryPosition*ones(1,nFrames);
end

for frame_index = 1:nFrames
    NumFurrowMeasurements(frame_index) = size( FurrowMeasurements{frame_index},1);
    %     for k =1: NumFurrowMeasurements(frame_index)
    %         FurrowMeasurements{frame_index}(k,3) = EmbryoBoundaryPositions(frame_index);
    %     end
    if NumFurrowMeasurements(frame_index) > 0
        
        EmbryoBoundaryPositions(frame_index) = mean(FurrowMeasurements{frame_index}(:,3));
        
    end
    
    NumPVMeasurements(frame_index) = size( PerivitellineMeasurements{frame_index},1);
    
end



%%
counter= 1;
CurrentFrame=nc14;
MemImage = membraneMat(:,:,CurrentFrame);
DisplayRange=[min(min(MemImage)),max(max(MemImage))];

Overlay=figure;
set(Overlay,'units', 'normalized', 'position',[0.01, .1, .8, .8]);

overlayAxes = axes(Overlay,'Units', 'normalized', 'Position', [0 0 1 1]);

tb = axtoolbar(overlayAxes);
tb.Visible = 'off';


currentCharacter=1;



% Show the first image
imOverlay = imshow(MemImage,DisplayRange,'Parent',overlayAxes);
hold(overlayAxes,'on')
x_ticks =[0:100:xSize xSize];
BoundaryPlotHandle = plot(overlayAxes, x_ticks, EmbryoBoundaryPositions(CurrentFrame)*ones(1,length(x_ticks)), 'b-+');
GuidePlotHandles = cell(1,length(x_ticks));
for pl=2:length(x_ticks)-1
    GuidePlotHandles{pl} = plot(overlayAxes, [x_ticks(pl),x_ticks(pl)], [0 ySize], 'b-+');
    GuidePlotHandles{pl}.Color(4) =0.3;
end
try
    FigureTitle=['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames),...
        ', nc: ',num2str(nc(CurrentFrame))];
catch
    FigureTitle=['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames)];
end
%%
title(overlayAxes, FigureTitle);


projFlag = false;
set(0, 'CurrentFigure', Overlay)


ShowFit =false;
PV = false;

%%

while (currentCharacter~='x')
    NumPVMeasurements(CurrentFrame) = size(PerivitellineMeasurements{CurrentFrame},1);
    NFurrowPoints=size(FurrowMeasurements{CurrentFrame}, 1);
    
    try
        FigureTitle=['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames),...
            ', nc: ',num2str(nc(CurrentFrame))];
    catch
        FigureTitle=['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames)];
    end
    
    set(Overlay,'Name',FigureTitle)
    MemImage = membraneMat(:,:,CurrentFrame);
    imOverlay.CData = MemImage;
    try
        caxis(overlayAxes, DisplayRange);
    end
    BoundaryPlotHandle.YData = EmbryoBoundaryPositions(CurrentFrame)*ones(1,length([0:100:xSize xSize]));
    
    
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
        
        PlotHandle{k} = plot(overlayAxes, FurrowMeasurements{CurrentFrame}(k,1),FurrowMeasurements{CurrentFrame}(k,2) ,...
            'r+', 'MarkerSize',40);
        LineHandle{k} = plot(overlayAxes,[FurrowMeasurements{CurrentFrame}(k,1),FurrowMeasurements{CurrentFrame}(k,1)],...
            [FurrowMeasurements{CurrentFrame}(k,2),FurrowMeasurements{CurrentFrame}(k,3)],...
            'r','LineStyle','-','Marker','+');
        
        
    end
    
    
    NPVPoints=size(PerivitellineMeasurements{CurrentFrame}, 1);
    if exist('PVPlotHandle', 'var')
        cellfun(@delete, PVPlotHandle);
    end
    if exist('PVLineHandle', 'var')
        cellfun(@delete, PVLineHandle);
    end
    PVPlotHandle = cell(NPVPoints, 1);
    PVLineHandle = cell(NPVPoints, 1);
    for k=1:NPVPoints
        
        PVPlotHandle{k} = plot(overlayAxes, PerivitellineMeasurements{CurrentFrame}(k,1),PerivitellineMeasurements{CurrentFrame}(k,2) ,...
            'y+', 'MarkerSize',40);
        PVLineHandle{k} = plot(overlayAxes,[PerivitellineMeasurements{CurrentFrame}(k,1),PerivitellineMeasurements{CurrentFrame}(k,1)],...
            [PerivitellineMeasurements{CurrentFrame}(k,2),PerivitellineMeasurements{CurrentFrame}(k,3)],...
            'y','LineStyle','-','Marker','+');
        
        
    end
    
    try
        FigureTitle=['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames),...
            ', nc: ',num2str(nc(CurrentFrame))];
    catch
        FigureTitle=['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames)];
    end
    
    
    %%
    
    tb = axtoolbar(overlayAxes);
    tb.Visible = 'off';
    
    
    ct=waitforbuttonpress;
    currentCharacter=get(Overlay,'currentcharacter');
    currentMouse=get(overlayAxes,'CurrentPoint');
    
    
    if (ct ~=0)&(currentCharacter=='s')
        for frame_index = 1:nFrames
            NumFurrowMeasurements(frame_index) = size( FurrowMeasurements{CurrentFrame},1);
            
            if NumFurrowMeasurements(frame_index) > 0
                DeltaFC_vector = abs(FurrowMeasurements{frame_index}(:,2).'-FurrowMeasurements{frame_index}(:,3).');
                DeltaFC_pixels(frame_index,1) = mean(DeltaFC_vector);
                DeltaFC_pixels(frame_index,2) = std(DeltaFC_vector);
                DeltaFC_um(frame_index,1) = mean(DeltaFC_vector*PixelSize_um);
                DeltaFC_um(frame_index,2) = std(DeltaFC_vector*PixelSize_um);
                EmbryoBoundaryPositions(frame_index) = mean(FurrowMeasurements{frame_index}(:,3));
            else
                DeltaFC_pixels(frame_index,1)= NaN;
                DeltaFC_pixels(frame_index,2)= NaN;
                DeltaFC_um(frame_index,1) = NaN;
                DeltaFC_um(frame_index,2) = NaN;
            end
            
            NumPVMeasurements(frame_index) = size( PerivitellineMeasurements{CurrentFrame},1);
            
            if NumPVMeasurements(frame_index) > 0
                DeltaPV_vector = abs(PerivitellineMeasurements{frame_index}(:,2).'-PerivitellineMeasurements{frame_index}(:,3).');
                DeltaPV_pixels(frame_index,1) = mean(DeltaPV_vector);
                DeltaPV_pixels(frame_index,2) = std(DeltaPV_vector);
                DeltaPV_um(frame_index,1) = mean(DeltaPV_vector*PixelSize_um);
                DeltaPV_um(frame_index,2) = std(DeltaPV_vector*PixelSize_um);
            else
                DeltaPV_pixels(frame_index,1)= NaN;
                DeltaPV_pixels(frame_index,2)= NaN;
                DeltaPV_um(frame_index,1) = NaN;
                DeltaPV_um(frame_index,2) = NaN;
            end
        end
        save([DropboxFolder,filesep,Prefix,filesep,'FurrowCanalDepthMeasurements.mat'],'FurrowMeasurements',...
            'NumFurrowMeasurements','DeltaFC_pixels','DeltaFC_um','EmbryoBoundaryPositions',...
            'PerivitellineMeasurements','NumPVMeasurements','DeltaPV_pixels', 'DeltaPV_um','-v6');
        
        
        
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'normal'))
        currentCharacter=1;
        if ~PV
        if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=ySize)&(currentMouse(1,1)<=xSize)
            
            %Add a circle to this location with the median radius of the
            %ellipses found in this frame
            
            %(x, y, a, b, theta, maxcontourvalue, time,
            
            [x,y] = ginputc(1, 'Color', 'r');
            FurrowMeasurements{CurrentFrame}(end+1,:)= [currentMouse(1,1),currentMouse(1,2), y];
            
            
        end
        else
        if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=ySize)&(currentMouse(1,1)<=xSize)
            
            %Add a circle to this location with the median radius of the
            %ellipses found in this frame
            
            %(x, y, a, b, theta, maxcontourvalue, time,
            
            [x,y] = ginputc(1, 'Color', 'y');
            PerivitellineMeasurements{CurrentFrame}(end+1,:)= [currentMouse(1,1),currentMouse(1,2), y];
            NumPVMeasurements(CurrentFrame) = size(PerivitellineMeasurements{CurrentFrame},1);
            
        end
        end
        
        
        
    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'alt'))
        currentCharacter=1;
        if ~PV
        if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=ySize)&(currentMouse(1,1)<=xSize)
            %Find out which ellipses we clicked on so we can delete it
            
            %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
            Distances=sqrt((FurrowMeasurements{CurrentFrame}(:,1)-currentMouse(1,1)).^2+...
                (FurrowMeasurements{CurrentFrame}(:,2)-currentMouse(1,2)).^2);
            if ~isempty(Distances)
                [~,MinIndex]=min(Distances);
                
                FurrowMeasurements{CurrentFrame}=[FurrowMeasurements{CurrentFrame}(1:MinIndex-1,:);...
                    FurrowMeasurements{CurrentFrame}(MinIndex+1:end,:)];
            end
        end
        else
            if (currentMouse(1,2)>0)&(currentMouse(1,1)>0)&(currentMouse(1,2)<=ySize)&(currentMouse(1,1)<=xSize)
            %Find out which ellipses we clicked on so we can delete it
            
            %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
            Distances=sqrt((PerivitellineMeasurements{CurrentFrame}(:,1)-currentMouse(1,1)).^2+...
                (PerivitellineMeasurements{CurrentFrame}(:,2)-currentMouse(1,2)).^2);
            if ~isempty(Distances)
                [~,MinIndex]=min(Distances);
                
                PerivitellineMeasurements{CurrentFrame}=[PerivitellineMeasurements{CurrentFrame}(1:MinIndex-1,:);...
                    PerivitellineMeasurements{CurrentFrame}(MinIndex+1:end,:)];
                NumPVMeasurements(CurrentFrame) = size(PerivitellineMeasurements{CurrentFrame},1);
            end
        end
        end
    
    elseif (ct~=0)&(currentCharacter=='.')&(CurrentFrame<nFrames)
        %          NumFurrowMeasurements(CurrentFrame) = size( FurrowMeasurements{CurrentFrame},1);
        %          if NumFurrowMeasurements(CurrentFrame) > 0
        %              for m=1:NumFurrowMeasurements(CurrentFrame)
        %                  FurrowMeasurements{CurrentFrame}(m,3) = EmbryoBoundaryPositions(CurrentFrame);
        %              end
        %          end
        
        CurrentFrame=CurrentFrame+1;
    elseif (ct~=0)&(currentCharacter=='>')&(CurrentFrame<nFrames)
        %          NumFurrowMeasurements(CurrentFrame) = size( FurrowMeasurements{CurrentFrame},1);
        %          if NumFurrowMeasurements(CurrentFrame) > 0
        %              for m=1:NumFurrowMeasurements(CurrentFrame)
        %                  FurrowMeasurements{CurrentFrame}(m,3) = EmbryoBoundaryPositions(CurrentFrame);
        %              end
        %          end
        if CurrentFrame + 4 < nFrames
        CurrentFrame=CurrentFrame+5;
        else
            CurrentFrame = nFrames;
        end
        
    elseif (ct~=0)&(currentCharacter==',')&(CurrentFrame>1)
        %           NumFurrowMeasurements(CurrentFrame) = size( FurrowMeasurements{CurrentFrame},1);
        %          if NumFurrowMeasurements(CurrentFrame) > 0
        %              for m=1:NumFurrowMeasurements(CurrentFrame)
        %                  FurrowMeasurements{CurrentFrame}(m,3) = EmbryoBoundaryPositions(CurrentFrame);
        %              end
        %          end
        
        CurrentFrame=CurrentFrame-1;
    elseif (ct~=0)&(currentCharacter=='<')&(CurrentFrame>nc14)
        %          NumFurrowMeasurements(CurrentFrame) = size( FurrowMeasurements{CurrentFrame},1);
        %          if NumFurrowMeasurements(CurrentFrame) > 0
        %              for m=1:NumFurrowMeasurements(CurrentFrame)
        %                  FurrowMeasurements{CurrentFrame}(m,3) = EmbryoBoundaryPositions(CurrentFrame);
        %              end
        %          end
        if CurrentFrame - 4 > nc14
        CurrentFrame=CurrentFrame-5;
        else
            CurrentFrame = nc14;
        end
    elseif (ct~=0)&(currentCharacter=='j')
        iJump=input('Frame to jump to: ');
        if (floor(iJump)>=nc14)&(iJump<=nFrames)
            CurrentFrame=iJump;
        else
            disp('Frame out of range.');
        end
        
    elseif (ct~=0)&(currentCharacter=='m')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)/1.5;
        
    elseif (ct~=0)&(currentCharacter=='n')    %Decrease contrast
        DisplayRange(2)=DisplayRange(2)*1.5;
        
    elseif (ct~=0)&(currentCharacter=='r')    %Reset the contrast
        DisplayRange=[min(min(MemImage)),max(max(MemImage))];
    elseif (ct~=0)&(currentCharacter=='f')    %Reset the contrast
        
        ShowFit = ~ShowFit;
    elseif (ct~=0)&(currentCharacter=='t')    %Switch Selction type between furrow and perivitelline space
        
        PV = ~PV;
    elseif (ct~=0)&(currentCharacter=='d')    %Delete all ellipses in the current frame
        if ~PV
            FurrowMeasurements{CurrentFrame}=NaN(0,3);
        else
           PerivitellineMeasurements{CurrentFrame}=NaN(0,3); 
        end
    elseif (ct~=0)&(currentCharacter == 'b')
        [x,y] = ginput;
        try
            p_b = polyfit(x.', y.',1);
            EmbryoBoundaryPositions(CurrentFrame)= mean(y);
        end
        %         NumFurrowMeasurements(CurrentFrame) = size( FurrowMeasurements{CurrentFrame},1);
        %         for k =1: NumFurrowMeasurements(CurrentFrame)
        %             FurrowMeasurements{CurrentFrame}(k,3) = EmbryoBoundaryPositions(CurrentFrame);
        %         end
    elseif (ct~=0)&(currentCharacter == 'c')&(CurrentFrame>nc14)
        if ~PV
        EmbryoBoundaryPositions(CurrentFrame) = EmbryoBoundaryPositions(CurrentFrame-1);
        else
            cp_index = find(NumPVMeasurements(1:CurrentFrame-1) > 0, 1, 'last');
            if ~isempty(cp_index)
            PerivitellineMeasurements{CurrentFrame} = PerivitellineMeasurements{cp_index};
            NumPVMeasurments(CurrentFrame) = size(PerivitellineMeasurements{CurrentFrame}, 1);
            end
        end
            
        %          NumFurrowMeasurements(CurrentFrame) = size( FurrowMeasurements{CurrentFrame},1);
        %         for k =1: NumFurrowMeasurements(CurrentFrame)
        %             FurrowMeasurements{CurrentFrame}(k,3) = EmbryoBoundaryPositions(CurrentFrame);
        %         end
    elseif (ct~=0)&(currentCharacter == 'v')&(CurrentFrame<nFrames)
        if ~PV
        EmbryoBoundaryPositions(CurrentFrame) = EmbryoBoundaryPositions(CurrentFrame+1);
        else
            cp_index = find(NumPVMeasurements(CurrentFrame+1:end) > 0, 1);
            if ~isempty(cp_index)
            PerivitellineMeasurements{CurrentFrame} = PerivitellineMeasurements{cp_index+CurrentFrame};
            NumPVMeasurments(CurrentFrame) = size(PerivitellineMeasurements{CurrentFrame}, 1);
            end
        end
    elseif (ct~=0)&(currentCharacter=='0')    %Debug mode
        keyboard
        
    end
end
%%
for frame_index = 1:nFrames
    NumFurrowMeasurements(frame_index) = size( FurrowMeasurements{frame_index},1);
    %     for k =1: NumFurrowMeasurements(frame_index)
    %         FurrowMeasurements{frame_index}(k,3) = EmbryoBoundaryPositions(frame_index);
    %     end
    if NumFurrowMeasurements(frame_index) > 0
        DeltaFC_vector = abs(FurrowMeasurements{frame_index}(:,2).'-FurrowMeasurements{frame_index}(:,3).');
        DeltaFC_pixels(frame_index,1) = mean(DeltaFC_vector);
        DeltaFC_pixels(frame_index,2) = std(DeltaFC_vector);
        DeltaFC_um(frame_index,1) = mean(DeltaFC_vector*PixelSize_um);
        DeltaFC_um(frame_index,2) = std(DeltaFC_vector*PixelSize_um);
        EmbryoBoundaryPositions(frame_index) = mean(FurrowMeasurements{frame_index}(:,3));
    else
        DeltaFC_pixels(frame_index,1)= NaN;
        DeltaFC_pixels(frame_index,2)= NaN;
        DeltaFC_um(frame_index,1) = NaN;
        DeltaFC_um(frame_index,2) = NaN;
    end
    NumPVMeasurements(frame_index) = size( PerivitellineMeasurements{frame_index},1);
            
            if NumPVMeasurements(frame_index) > 0
                DeltaPV_vector = abs(PerivitellineMeasurements{frame_index}(:,2).'-PerivitellineMeasurements{frame_index}(:,3).');
                DeltaPV_pixels(frame_index,1) = mean(DeltaPV_vector);
                DeltaPV_pixels(frame_index,2) = std(DeltaPV_vector);
                DeltaPV_um(frame_index,1) = mean(DeltaPV_vector*PixelSize_um);
                DeltaPV_um(frame_index,2) = std(DeltaPV_vector*PixelSize_um);
            else
                DeltaPV_pixels(frame_index,1)= NaN;
                DeltaPV_pixels(frame_index,2)= NaN;
                DeltaPV_um(frame_index,1) = NaN;
                DeltaPV_um(frame_index,2) = NaN;
            end
    
end
EarliestFrames = true;
for frame_index = nc14:nFrames
    if NumPVMeasurements(frame_index) == 0 & EarliestFrames
        pv_index = find(NumPVMeasurements(frame_index+1:end) > 0, 1);
        if ~isempty(pv_index)
        PerivitellineMeasurements{frame_index} = PerivitellineMeasurements{pv_index+frame_index};
        NumPVMeasurements(frame_index) = size(PerivitellineMeasurements{frame_index}, 1);
        DeltaPV_vector = abs(PerivitellineMeasurements{frame_index}(:,2).'-PerivitellineMeasurements{frame_index}(:,3).');
                DeltaPV_pixels(frame_index,1) = mean(DeltaPV_vector);
                DeltaPV_pixels(frame_index,2) = std(DeltaPV_vector);
                DeltaPV_um(frame_index,1) = mean(DeltaPV_vector*PixelSize_um);
                DeltaPV_um(frame_index,2) = std(DeltaPV_vector*PixelSize_um);
        end
    elseif NumPVMeasurements(frame_index) == 0 & ~EarliestFrames
        pv_index = find(NumPVMeasurements(frame_index-5:frame_index-1) > 0, 1, 'last');
        if ~isempty(pv_index)
        PerivitellineMeasurements{frame_index} = PerivitellineMeasurements{pv_index+frame_index-5-1};
        PerivitellineMeasurements{frame_index} = PerivitellineMeasurements{pv_index+frame_index-5-1};
        DeltaPV_vector = abs(PerivitellineMeasurements{frame_index}(:,2).'-PerivitellineMeasurements{frame_index}(:,3).');
                DeltaPV_pixels(frame_index,1) = mean(DeltaPV_vector);
                DeltaPV_pixels(frame_index,2) = std(DeltaPV_vector);
                DeltaPV_um(frame_index,1) = mean(DeltaPV_vector*PixelSize_um);
                DeltaPV_um(frame_index,2) = std(DeltaPV_vector*PixelSize_um);
        end
    else
        EarliestFrames = false;
    end
    
    
end
DeltaFCnoPV_um = NaN(size(DeltaFC_um));
DeltaFCnoPV_um(:,1) = DeltaFC_um(:,1)-DeltaPV_um(:,1);
DeltaFCnoPV_um(:,2) = sqrt(DeltaFC_um(:,2).^2+DeltaPV_um(:,2).^2);
DeltaFCnoPV_pixels = NaN(size(DeltaFC_pixels));
DeltaFCnoPV_pixels(:,1) = DeltaFC_pixels(:,1)-DeltaPV_pixels(:,1);
DeltaFCnoPV_pixels(:,2) = sqrt(DeltaFC_pixels(:,2).^2+DeltaPV_pixels(:,2).^2);
save([DropboxFolder,filesep,Prefix,filesep,'FurrowCanalDepthMeasurements.mat'],'FurrowMeasurements',...
            'NumFurrowMeasurements','DeltaFC_pixels','DeltaFC_um','EmbryoBoundaryPositions',...
            'PerivitellineMeasurements','NumPVMeasurements','DeltaPV_pixels', 'DeltaPV_um',...
            'DeltaFCnoPV_um','DeltaFCnoPV_pixels','-v6');
        
        
close all
figure(2)
plot(FrameTimes_NC14, DeltaFC_um(nc14:end,1).','k.')
hold on 
plot(FrameTimes_NC14, DeltaFCnoPV_um(nc14:end,1).','r.')
%plot(0:length(DeltaFC_um(nc14:end,1))-1, DeltaFC_um(nc14:end,1).','k.')
xlabel('Time into NC14 (min)')
ylabel('Depth of Furrow Canal (um)')
xlim([0 max(FrameTimes_NC14)])
title('Membrane Furrow Plot')
