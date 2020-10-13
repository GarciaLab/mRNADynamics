function [coordT, coordB, coordL, coordR] = SelectCytoplasmRegion(imm)
close all
coordL = 1;
coordR = size(imm, 2);
coordT = 1;
coordB = size(imm, 1);
DisplayRange=[min(min(imm)),max(max(imm))];
HistoneImage=figure;
apAx = axes(HistoneImage);
%Now, do the correction
cc=1;

while (cc~='x')
    
    
    imshow(imm,DisplayRange, 'Parent', apAx)
    %imshow(APImage,DisplayRange)
    
    axis image
    %axis off
    title({'left boundary (green), right boundary (red),', 'top boundary (blue), bottom boundary (yellow); original'})
    hold on
    
    if exist('coordL', 'var')
        xline(coordL,'g')
    end
    
    if exist('coordR', 'var')
        xline(coordR,'r')
    end
    
    if exist('coordT', 'var')
        yline(coordT,'b')
    end
    
    if exist('coordB', 'var')
        yline(coordB,'y')
    end
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
    
    figure(HistoneImage)
    ct=waitforbuttonpress;
    cc=get(HistoneImage,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    
    if (ct~=0)&(cc=='c')        %Clear all AP information
        leftboundary=[];
        rightboundary=[];
        topboundary=[];
        bottomboundary=[];
    elseif (ct~=0)&(cc=='l')	%Select anterior end
        [coordL,coordLy]=ginputc(1,'Color',[1,1,1]);
    elseif (ct~=0)&(cc=='r')    %Select posterior end
        [coordR,coordRy]=ginputc(1,'Color',[1,1,1]);
    elseif (ct~=0)&(cc=='t')    %Select posterior end
        [coordTx,coordT]=ginputc(1,'Color',[1,1,1]);
    elseif (ct~=0)&(cc=='b')    %Select posterior end
        [coordBx,coordB]=ginputc(1,'Color',[1,1,1]);
    end
end

close all

coordB = uint16(coordB);
coordT = uint16(coordT);
coordR = uint16(coordR);
coordL = uint16(coordL);
end
