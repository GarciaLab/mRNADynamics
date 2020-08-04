function [ShiftColumn,ShiftRow]=ManualAPCorrection(SurfImage,ZoomImage,C,ResizeFactor,ShiftRow,ShiftColumn,  FullEmbryo, ZoomRatio, SurfRows, Rows, Columns, coordA, coordP, SurfColumns)
%
% RowsZoom/2+*2+RowsZoom/2-1+,...
%         2-ColumnsZoom/2+ShiftColumn:/2+ColumnsZoom/);
%
%
% TopLeftHalf,BottomRightHalf,...
%        ShiftColumn,ShiftRow,coordAHalf,coordPHalf,Rows,Columns,Zoom)
%
%m - Place the area with your mouse
%. - Move to the right
%> - Move to the right further
%, - Move to the left
%< - Move to the left further
%a - Move up
%A - Move up further
%z - Move down
%Z - Move down further
%x - Save and cancel

%Make an overlay of the zoomed in and zoomed out real
%images as well as of a quickly segmented nuclear mask

%Close existing images
close all

%Information about the correlation image
[CRows,CColumns]=size(C);

%Resize the zoom out
ZoomOutResized=imresize(SurfImage, ResizeFactor);

%Get the nuclear mask overlay
NucMaskZoomOut=GetNuclearMask(SurfImage,2.5,0);
NucMaskZoomOutResized=imresize(NucMaskZoomOut, ResizeFactor);
NucMaskZoomIn=GetNuclearMask(ZoomImage,8,2);

%Sizes of images
[RowsResized,ColumnsResized]=size(NucMaskZoomOutResized);
[RowsZoom,ColumnsZoom]=size(ZoomImage);

FigureCorrelation=figure;
axesCorrelation = axes(FigureCorrelation);
set(FigureCorrelation,'units', 'normalized', 'position',[0.4, 0.2, 0.3, 0.2]);
contourf(imresize(abs(C),.1)); %AR 9/12/2018. resized the image to make this more computationally realistic and prevent crashes. 
ylim(axesCorrelation,[(CRows-1)/2-RowsZoom,(CRows-1)/2+RowsZoom])
xlim(axesCorrelation,[(CColumns-1)/2-ColumnsZoom,(CColumns-1)/2+ColumnsZoom])
PlotHandle=[];

FigureRealOverlay=figure;
realAxes = axes(FigureRealOverlay);
set(FigureRealOverlay,'units', 'normalized', 'position',[0.05, 0.2, 0.3, 0.3]);

FigureOverlay=figure;
overlayAxes = axes(FigureOverlay);
set(FigureOverlay,'units', 'normalized', 'position',[0.05, 0.6, 0.3, 0.3]);

cc=1;

FigureRectangle=figure;
rectAx = axes(FigureRectangle);
set(FigureRectangle,'units', 'normalized', 'position',[0.4, 0.5, 0.3, 0.4]);


%Default flags
Green=true;
Red=true;

%Overlay the zoom in and zoom out images
while (cc~=13)
    

    %Show the correlation image and the shift position we're at
    delete(PlotHandle)
    hold(axesCorrelation, 'on')
    CorrX=ShiftColumn*ResizeFactor+(CColumns/2+1);
    CorrY=ShiftRow*ResizeFactor+(CRows/2+1);
    PlotHandle=plot(axesCorrelation,CorrX,CorrY,...
        'or','MarkerSize',10);
    hold(axesCorrelation, 'off')

   
    
    %Crop the zoomed out nuclear mask
    rowInd1 = uint16( (RowsResized/2-RowsZoom/2) + (ShiftRow*ResizeFactor) + 1 );
    rowInd2 = uint16(round(RowsResized/2+RowsZoom/2+ShiftRow*ResizeFactor));
    colInd1 = uint16( (ColumnsResized/2-ColumnsZoom/2) + (ShiftColumn*ResizeFactor) + 1);
    colInd2 = uint16(ColumnsResized/2+ColumnsZoom/2+ShiftColumn*ResizeFactor);
    try
        NucMaskZoomOutResizedCropped = NucMaskZoomOutResized(rowInd1:rowInd2, colInd1:colInd2);
    
        ZoomOutResizedCropped=...
            ZoomOutResized(rowInd1:rowInd2,colInd1:colInd2);
    
   
        ImOverlayMask=cat(3,mat2gray(NucMaskZoomOutResizedCropped),...
            +mat2gray(NucMaskZoomIn),zeros(size(NucMaskZoomOutResizedCropped)));
        ImOverlay=cat(3,mat2gray(ZoomOutResizedCropped)*Red,...
            +mat2gray(ZoomImage)*Green,zeros(size(ZoomOutResizedCropped)));

        imshow(ImOverlayMask,'Parent', overlayAxes)
        imshow(ImOverlay, 'Parent', realAxes);
    catch
        warning('Could not generate correlation image.');
    end
    set(FigureRectangle,'name',['ShiftRow: ',num2str(ShiftRow),'. ShiftColumn:',num2str(ShiftColumn),'.'])


    ImageCenter=[SurfRows/2 + ShiftRow,SurfColumns/2 + ShiftColumn];
    TopLeft=[ImageCenter(1)-Rows/ZoomRatio/2,ImageCenter(2)-Columns/ZoomRatio/2];
    BottomRight=[ImageCenter(1)+Rows/ZoomRatio/2,ImageCenter(2)+Columns/ZoomRatio/2];
    imshow(imadjust(mat2gray(SurfImage)),'DisplayRange',[],'InitialMagnification',100, 'Parent', rectAx)
    hold(rectAx, 'on')
    rectangle(rectAx,'Position',[TopLeft([2,1]),BottomRight([2,1])-TopLeft([2,1])],'EdgeColor','r')
    plot(rectAx,coordA(1),coordA(2),'.g','MarkerSize',30)
    plot(rectAx,coordP(1),coordP(2),'.r','MarkerSize',30)
    plot(rectAx,[coordA(1),coordP(1)],[coordA(2),coordP(2)],'-b')
    hold(rectAx, 'off')

    
%     figure(FigureOverlay)
    figure(FigureRectangle)
    drawnow
    ct = waitforbuttonpress;
    ct=1;
    cc=get(FigureRectangle,'currentcharacter');
    %cm=get(gca,'CurrentPoint');
    

    %Move up/down + left/right
    if (ct~=0)&(cc=='.')
        ShiftColumn=ShiftColumn+1;
    elseif (ct~=0)&(cc==',')
        ShiftColumn=ShiftColumn-1;
    elseif (ct~=0)&(cc=='z')
        ShiftRow=ShiftRow+1;
    elseif (ct~=0)&(cc=='a')
        ShiftRow=ShiftRow-1;
    elseif (ct~=0)&(cc=='>')
        ShiftColumn=ShiftColumn+50;
    elseif (ct~=0)&(cc=='<')
        ShiftColumn=ShiftColumn-50;
    elseif (ct~=0)&(cc=='Z')
        ShiftRow=ShiftRow+50;
    elseif (ct~=0)&(cc=='A')
        ShiftRow=ShiftRow-50;
    elseif (ct~=0)&(cc=='m')
        disp('Select a location to center the alignment box');
        [Positionx,Positiony]=ginputc(1,'color', 'b', 'linewidth',1);
        Position = [Positionx,Positiony];
        ShiftColumn=round(Position(1)-SurfColumns/2);
        ShiftRow=round(Position(2)-SurfRows/2);
        
       
    %Turn channels on and off
    elseif (ct~=0)&(cc=='r')
        Red=~Red;
    elseif (ct~=0)&(cc=='g')
        Green=~Green;
        
    
    %Zoom in and out
    elseif (ct~=0)&(cc=='=') %#ok<*AND2>
        OldXLim=xlim(axesCorrelation);
        OldYLim=ylim(axesCorrelation);
        xlim(axesCorrelation,[OldXLim-CorrX]/2+CorrX)
        ylim(axesCorrelation,[OldYLim-CorrY]/2+CorrY)
       
        
    elseif (ct~=0)&(cc=='-')
        OldXLim=xlim(axesCorrelation);
        OldYLim=ylim(axesCorrelation);
        xlim(axesCoxrrelation,[OldXLim-CorrX]*2+CorrX)
        ylim(axesCorrelation,[OldYLim-CorrY]*2+CorrY)
       
        
    elseif (ct~=0)&(cc=='x')
        cc=13;
    elseif (ct~=0)&(cc=='9')
        keyboard
    end
end

close(FigureOverlay);
close(FigureRectangle);
close(FigureRealOverlay);
close(FigureCorrelation);