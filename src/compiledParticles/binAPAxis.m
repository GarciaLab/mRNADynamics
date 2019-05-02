function [APbinID, APbinArea] = binAPAxis(APResolution, ...
    FrameInfo, coordAZoom, APAngle, APLength, minBinSize)
%binAPAxis Summary of this function goes here
%   Detailed explanation goes here

%minBinSize is the fraction of the median bin size;

APbinID=0:APResolution:1;

%Create an image for the different AP bins
APPosImage=zeros(FrameInfo(1).LinesPerFrame,FrameInfo(1).PixelsPerLine);
[Rows,Columns]=size(APPosImage);

for i=1:Rows
    for j=1:Columns
        try
            Angle=atan((i-coordAZoom(2))./(j-coordAZoom(1)));
        catch
            Angle=atan2((i-coordAZoom(2)),(j-coordAZoom(1)));
        end
        Distance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
        APPosition=Distance.*cos(Angle-APAngle);
        APPosImage(i,j)=APPosition/APLength;
    end
end


APPosBinImage=zeros(size(APPosImage));
for i=1:(length(APbinID)-1)
    FilteredMask=(APbinID(i)<=APPosImage)&(APbinID(i+1)>APPosImage);
    APPosBinImage=APPosBinImage+FilteredMask*i;
end


%Calculate the area in pixels corresponding to each AP bin. We will use
%this to get rid of small AP bins in the image and also to calculate
%probabilities of nuclei being active.
APbinArea = zeros(length(APbinID),1);
%Calculate ther areas of the AP bins
for i=1:length(APbinID)
    APbinArea(i)=sum(sum(APPosBinImage==i));
end
%Get the median of the non-zero areas
MedianArea=median(APbinArea(APbinArea>0));
%Only keep the bins with an area of at least minBinSize of the median
APbinArea(APbinArea<MedianArea*minBinSize)= NaN;
end

