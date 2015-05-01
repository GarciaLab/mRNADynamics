function [ShiftColumn,ShiftRow]=ManualAPCorrection(SurfImage,ZoomImage,C,ResizeFactor,ShiftRow,ShiftColumn, FullEmbryo, ZoomRatio, SurfRows, Rows, Columns, coordA, coordP, SurfColumns)

%. - Move to the right
%, - Move to the left
%a - Move up
%z - Move down

FigureRealOverlay=figure;
set(gcf,'Position',[22  -212   682   347])

FigureOverlay=figure;
set(gcf,'Position',[15   324   676   342])
cc=1;

%Default flags
Green=true;
Red=true;

while (cc~=13)
    
    %%%
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
    
    %%%
    figure(FigureOverlay)
    ct=waitforbuttonpress;
    cc=get(FigureOverlay,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    set(gcf,'name',['ShiftRow: ',num2str(ShiftRow),'. ShiftColumn:',num2str(ShiftColumn),'.'])

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
        ShiftColumn=ShiftColumn+30;
    elseif (ct~=0)&(cc=='<')
        ShiftColumn=ShiftColumn-30;
    elseif (ct~=0)&(cc=='A')
        ShiftRow=ShiftRow-30;
    elseif (ct~=0)&(cc=='Z')
        ShiftRow=ShiftRow+30;
        
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



