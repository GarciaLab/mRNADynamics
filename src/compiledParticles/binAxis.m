function [binID, binArea] = binAxis(axisResolution, ...
    FrameInfo, coordAZoom, axisAngle, axisLength, minBinSize, axis)

if strcmpi(axis, 'AP')
    binID=0:axisResolution:1;
elseif strcmpi(axis, 'DV')
    binID = 0:axisResolution:1000; %why -800? AR 6/3/2019
end

%Create an image for the different AP bins
posImage=zeros(FrameInfo(1).LinesPerFrame,FrameInfo(1).PixelsPerLine);
[Rows,Columns]=size(posImage);

for i=1:Rows
    for j=1:Columns
        try
            Angle=atan((i-coordAZoom(2))./(j-coordAZoom(1)));
        catch
            Angle=atan2((i-coordAZoom(2)),(j-coordAZoom(1)));
        end
        Distance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
        if strcmpi(axis, 'AP')
            position=Distance.*cos(Angle-axisAngle);
            posImage(i,j)=position/axisLength;
        elseif strcmpi(axis, 'DV') 
            position = Distance.*sin(Angle-axisAngle);
            posImage(i,j) = position;
        end
    end
end


posBinImage=zeros(size(posImage));
for i=1:(length(binID)-1)
    FilteredMask=(binID(i)<=posImage)&(binID(i+1)>posImage);
    posBinImage=posBinImage+FilteredMask*i;
end


%Calculate the area in pixels corresponding to each bin. We will use
%this to get rid of small bins in the image and also to calculate
%probabilities of nuclei being active.
binArea = zeros(length(binID),1);
%Calculate ther areas of the AP bins
for i=1:length(binID)
    binArea(i)=sum(sum(posBinImage==i));
end
%Get the median of the non-zero areas
MedianArea=median(binArea(binArea>0));
%Only keep the bins with an area of at least minBinSize of the median
binArea(binArea<MedianArea*minBinSize)= NaN;


end