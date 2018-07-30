function [x,y,z]=SpotsXYZ(Spots)

%Return the X and Y coordinate of the brightest Z of each spot in the
%Spots structure

if ~isempty(Spots.Fits)
    for i=1:length(Spots.Fits)
        x(i)=Spots.Fits(i).xDoG(Spots.Fits(i).z==Spots.Fits(i).brightestZ);
        y(i)=Spots.Fits(i).yDoG(Spots.Fits(i).z==Spots.Fits(i).brightestZ);
        z(i)=Spots.Fits(i).brightestZ;
    end
else
    x=[];
    y=[];
    z=[];
end

