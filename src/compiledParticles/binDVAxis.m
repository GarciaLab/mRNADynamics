function [DVbinID, DVbinArea] = binDVAxis(FrameInfo, coordAZoom, APAngle)
%binDVAxis Summary of this function goes here
%   Detailed explanation goes here

DVbinID=linspace(-800,0,51); %JAKE: Would change to DV resolution later
%Create an image for the different DV bins
DVPosImage=zeros(FrameInfo(1).LinesPerFrame,FrameInfo(1).PixelsPerLine);
[Rows,Columns]=size(DVPosImage);

for i=1:Rows
    for j=1:Columns
        Angle=atan2((i-coordAZoom(2)),(j-coordAZoom(1)));   

        Distance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
        DVPosition=Distance.*sin(Angle-APAngle);
        DVPosImage(i,j)=DVPosition;
    end
end


DVPosBinImage=zeros(size(DVPosImage));
for i=1:(length(DVbinID)-1)
    FilteredMask=(DVbinID(i)<=DVPosImage)&(DVbinID(i+1)>DVPosImage);
    DVPosBinImage=DVPosBinImage+FilteredMask*i;
end


%Calculate the area in pixels corresponding to each AP bin. We will use
%this to get rid of small AP bins in the image and also to calculate
%probabilities of nuclei being active.
DVbinArea = zeros(length(DVbinID),1);
%Calculate ther areas of the AP bins
for i=1:length(DVbinID)
    DVbinArea(i)=sum(sum(DVPosBinImage==i));
end
%Get the median of the non-zero areas
MedianArea=median(DVbinArea(DVbinArea>0));
%Only keep the bins with an area of at least 70% of the median
DVbinArea(DVbinArea<MedianArea*0.7)=nan;
end

