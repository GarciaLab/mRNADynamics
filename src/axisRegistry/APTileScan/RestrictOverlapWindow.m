% RestrictOverlapWindow.m
% author: Gabriella Martini
% date created: 7/30/20
% date last modified: 7/30/20



APImageFig=figure;
apAx = axes(APImageFig);
%Now, do the correction
cc=1;

while (cc~='x')
    
    
    imshow(imagesc(C), 'Parent', apAx)
    %imshow(APImage,DisplayRange)
    
    axis image
    axis off
    title('left boundary (green), right boundary (red), top boundary (blue), bottom boundary (yellow); original')
    hold on
%     
%     try
%         plot(coordA(1),coordA(2),'g.','MarkerSize',20);
%     catch
%         %not sure what happened here.
%     end
%     
%     try
%         plot(coordP(1),coordP(2),'r.','MarkerSize',20);
%     catch
%         %not sure what happened here.
%     end
    

    
    hold off
    
    figure(APImageFig)
    ct=waitforbuttonpress;
    cc=get(APImageFig,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    
    if (ct~=0)&(cc=='c')        %Clear all AP information
        leftboundary=[];
        rightboundary=[];
        topboundary=[];
        bottomboundary=[];
    elseif (ct~=0)&(cc=='a')	%Select anterior end
        [coordAx,CoordAy]=ginputc(1,'Color',[1,1,1]);
        coordA = [coordAx,CoordAy];
    elseif (ct~=0)&(cc=='p')    %Select posterior end
        [coordPx,CoordPy]=ginputc(1,'Color',[1,1,1]);
        coordP = [coordPx,CoordPy];
    elseif (ct~=0)&(cc=='.')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)/2;
        
    elseif (ct~=0)&(cc==',')    %Decrease contrast
        DisplayRange(2)=DisplayRange(2)*2;
        
    elseif (ct~=0)&(cc=='r')    %Reset the contrast
        DisplayRange=[min(min(APImage)),max(max(APImage))];
        
    elseif (ct==0)&(strcmp(get(APImageFig,'SelectionType'),'alt')) %Delete the point that was clicked on
        cc=1;
        
        [~,MinIndex]=min((cm(1,1)-[coordA(1),coordP(1)]).^2+(cm(1,2)-[coordA(2),coordP(2)]).^2);
        
        if MinIndex==1
            coordA=[];
        elseif MinIndex==2
            coordP=[];
        end
    elseif (ct~=0)&(cc=='m')        %Manual stitching mode
        %ManualStitch
    elseif (ct~=0)&(cc=='s')
        coordPTemp=coordA;
        coordA=coordP;
        coordP=coordPTemp;
    elseif (ct~=0)&(cc=='v')
        [coordVx,CoordVy]=ginputc(1,'Color','m');
        coordV = [coordVx,CoordVy];
        saveVars = [saveVars, 'coordV'];
        dv = true;
    elseif (ct~=0)&(cc=='d')
        [coordDx,CoordDy]=ginputc(1,'Color','y');
        coordD = [coordDx,CoordDy];
        saveVars = [saveVars, 'coordD'];
        dv = true;
        
        
    end
end
