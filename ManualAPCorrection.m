function [ShiftColumn,ShiftRow,TopLeftHalf,BottomRightHalf]=ManualAPCorrection(SurfImage,TopLeftHalf,BottomRightHalf,...
       ShiftColumn,ShiftRow,coordAHalf,coordPHalf,Rows,Columns,Zoom)

%. - Move to the right
%, - Move to the left
%a - Move up
%z - Move down


FigureOverlay=figure;
cc=1;

%Plot the image and overlay with the different particles found

while (cc~=13)
    
    [NRows,NCols]=size(SurfImage);
    ImageCenter=[NRows/2,NCols/2];
    
    TopLeftHalf=[ImageCenter(1)-Rows/Zoom/2+ShiftRow,...
        ImageCenter(2)-Columns/Zoom/2+ShiftColumn];
    BottomRightHalf=[ImageCenter(1)+Rows/Zoom/2+ShiftRow,...
        ImageCenter(2)+Columns/Zoom/2+ShiftColumn];

    
    imshow(SurfImage,'DisplayRange',[],'InitialMagnification',100)
    hold on
    rectangle('Position',[TopLeftHalf([2,1]),BottomRightHalf([2,1])-TopLeftHalf([2,1])],'EdgeColor','r')
    plot(coordAHalf(1),coordAHalf(2),'.g','MarkerSize',30)
    plot(coordPHalf(1),coordPHalf(2),'.r','MarkerSize',30)
    plot([coordAHalf(1),coordPHalf(1)],[coordAHalf(2),coordPHalf(2)],'-b')
    hold off
    xlim([NCols/2*0.6,NCols/2*1.4])
    ylim([NRows/2*0.6,NRows/2*1.4])
    
    
    
    set(gcf,'name',['ShiftRow: ',num2str(ShiftRow),'. ShiftColumn:',num2str(ShiftColumn),'.'])
    
    figure(FigureOverlay)
    ct=waitforbuttonpress;
    cc=get(FigureOverlay,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    


    if (ct~=0)&(cc=='.')
        ShiftColumn=ShiftColumn+1;
    elseif (ct~=0)&(cc==',')
        ShiftColumn=ShiftColumn-1;
    elseif (ct~=0)&(cc=='a')
        ShiftRow=ShiftRow-1;
    elseif (ct~=0)&(cc=='z')
        ShiftRow=ShiftRow+1;
    end
end

close(FigureOverlay)
