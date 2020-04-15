function NewImage=ShiftImage(Image,ShiftX,ShiftY)

%Shifts Image by ShiftX and ShiftY. It keep the size of the image constant
%by padding with zeros.

[Rows,Columns]=size(Image);

NewImage=zeros(size(Image));
%NewImage=ones(size(Image))*0.5;

%This is because I screwed up somewhere and I'm too lazy to fix it right
%now.
ShiftX=-ShiftX;
ShiftY=-ShiftY;



if ShiftX<0
    OriginalRangeX=(-ShiftX+1):Columns;
    TargetRangeX=1:(Columns+ShiftX);
elseif ShiftX>0
    OriginalRangeX=1:(Columns-ShiftX);
    TargetRangeX=(ShiftX+1):(Columns);
else
    OriginalRangeX=1:Columns;
    TargetRangeX=1:Columns;
end


if ShiftY<0
    OriginalRangeY=(-ShiftY+1):Rows;
    TargetRangeY=1:(Rows+ShiftY);
elseif ShiftY>0
    OriginalRangeY=1:(Rows-ShiftY);
    TargetRangeY=(ShiftY+1):(Rows);
else
    OriginalRangeY=1:Rows;
    TargetRangeY=1:Rows;
end




NewImage(TargetRangeY,TargetRangeX)=Image(OriginalRangeY,OriginalRangeX);
% 
% figure(1)
% imshow(Image)
% figure(2)
% imshow(NewImage)
% 
