function [ShiftColumn,ShiftRow]=ManualAPCorrection(SurfImage,ZoomImage,C,ResizeFactor,ShiftRow,ShiftColumn,  FullEmbryo, ZoomRatio, SurfRows, Rows, Columns, coordA, coordP, SurfColumns)

% RowsZoom/2+*2+RowsZoom/2-1+,...
%         2-ColumnsZoom/2+ShiftColumn:/2+ColumnsZoom/);


% TopLeftHalf,BottomRightHalf,...
%        ShiftColumn,ShiftRow,coordAHalf,coordPHalf,Rows,Columns,Zoom)

%. - Move to the right
%, - Move to the left
%a - Move up
%z - Move down
%x - Save and cancel


%Make an overlay of the zoomed in and zoomed out real
%images as well as of a quickly segmented nuclear mask

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
set(gcf,'Position',[709   263   560   420])
contourf(abs(C))
ylim([(CRows-1)/2-RowsZoom,(CRows-1)/2+RowsZoom])
xlim([(CColumns-1)/2-ColumnsZoom,(CColumns-1)/2+ColumnsZoom])
PlotHandle=[];

FigureRealOverlay=figure;
set(gcf,'Position',[22  -212   682   347])

FigureOverlay=figure;
set(gcf,'Position',[15   324   676   342])
cc=1;

FigureRectangle=figure;
% set(gcf,'Position',[115   424   776   442])


%Default flags
Green=true;
Red=true;

%Overlay the zoom in and zoom out images
while (cc~=13)
    

    %Show the correlation image and the shift position we're at
    figure(FigureCorrelation)
    delete(PlotHandle)
    hold on
    CorrX=ShiftColumn*ResizeFactor+(CColumns/2+1);
    CorrY=ShiftRow*ResizeFactor+(CRows/2+1);
    PlotHandle=plot(CorrX,CorrY,...
        'or','MarkerSize',10);
    hold off
  
   
    
    %Crop the zoomed out nuclear mask
    NucMaskZoomOutResizedCropped=...
        NucMaskZoomOutResized(RowsResized/2-RowsZoom/2+ShiftRow*ResizeFactor+1:RowsResized/2+RowsZoom/2+ShiftRow*ResizeFactor,...
        ColumnsResized/2-ColumnsZoom/2+ShiftColumn*ResizeFactor+1:ColumnsResized/2+ColumnsZoom/2+ShiftColumn*ResizeFactor);
    ZoomOutResizedCropped=...
        ZoomOutResized(RowsResized/2-RowsZoom/2+ShiftRow*ResizeFactor+1:RowsResized/2+RowsZoom/2+ShiftRow*ResizeFactor,...
        ColumnsResized/2-ColumnsZoom/2+ShiftColumn*ResizeFactor+1:ColumnsResized/2+ColumnsZoom/2+ShiftColumn*ResizeFactor);
    
   
    ImOverlayMask=cat(3,mat2gray(NucMaskZoomOutResizedCropped),...
        +mat2gray(NucMaskZoomIn),zeros(size(NucMaskZoomOutResizedCropped)));
    ImOverlay=cat(3,mat2gray(ZoomOutResizedCropped)*Red,...
        +mat2gray(ZoomImage)*Green,zeros(size(ZoomOutResizedCropped)));
    
    figure(FigureOverlay)
    imshow(ImOverlayMask)
    set(gcf,'name',['ShiftRow: ',num2str(ShiftRow),'. ShiftColumn:',num2str(ShiftColumn),'.'])

    
    figure(FigureRealOverlay)
    imshow(ImOverlay)
    
    figure(FigureRectangle)
    ImageCenter=[SurfRows/2 + ShiftRow,SurfColumns/2 + ShiftColumn];
    TopLeft=[ImageCenter(1)-Rows/ZoomRatio/2,ImageCenter(2)-Columns/ZoomRatio/2];
    BottomRight=[ImageCenter(1)+Rows/ZoomRatio/2,ImageCenter(2)+Columns/ZoomRatio/2];
    imshow(imadjust(mat2gray(FullEmbryo)),'DisplayRange',[],'InitialMagnification',100)
    hold on
    rectangle('Position',[TopLeft([2,1]),BottomRight([2,1])-TopLeft([2,1])],'EdgeColor','r')
    plot(coordA(1),coordA(2),'.g','MarkerSize',30)
    plot(coordP(1),coordP(2),'.r','MarkerSize',30)
    plot([coordA(1),coordP(1)],[coordA(2),coordP(2)],'-b')
    hold off
    
    figure(FigureOverlay)
    ct=waitforbuttonpress;
    cc=get(FigureOverlay,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    

    %Move up/down + left/right
    if (ct~=0)&(cc=='.')
        ShiftColumn=ShiftColumn+1;
    elseif (ct~=0)&(cc==',')
        ShiftColumn=ShiftColumn-1;
    elseif (ct~=0)&(cc=='a')
        ShiftRow=ShiftRow-1;
    elseif (ct~=0)&(cc=='z')
        ShiftRow=ShiftRow+1;
    elseif (ct~=0)&(cc=='>')
        ShiftColumn=ShiftColumn+10;
    elseif (ct~=0)&(cc=='<')
        ShiftColumn=ShiftColumn-10;
    elseif (ct~=0)&(cc=='A')
        ShiftRow=ShiftRow-10;
    elseif (ct~=0)&(cc=='Z')
        ShiftRow=ShiftRow+10;
        
    %Turn channels on and off
    elseif (ct~=0)&(cc=='r')
        Red=~Red;
    elseif (ct~=0)&(cc=='g')
        Green=~Green;
        
    
    %Zoom in and out
    elseif (ct~=0)&(cc=='=')
        figure(FigureCorrelation)
        OldXLim=xlim;
        OldYLim=ylim;
        xlim([OldXLim-CorrX]/2+CorrX)
        ylim([OldYLim-CorrY]/2+CorrY)
    elseif (ct~=0)&(cc=='-')
        figure(FigureCorrelation)
        OldXLim=xlim;
        OldYLim=ylim;
        xlim([OldXLim-CorrX]*2+CorrX)
        ylim([OldYLim-CorrY]*2+CorrY)
    elseif (ct~=0)&(cc=='x')
        cc=13;
    elseif (ct~=0)&(cc=='9')
        keyboard
    end
end

close(FigureOverlay)
close(FigureRectangle)
close(FigureRealOverlay)



