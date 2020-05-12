function [x,y,z]=SpotsXYZ(Spots, varargin)

useGauss3DCentroid = false;
if ~isempty(varargin)
    useGauss3DCentroid = true;
end
%Return the X and Y coordinate of the brightest Z of each spot in the
%Spots structure

if ~isempty(Spots.Fits)
    for f=1:length(Spots.Fits)
        if ~useGauss3DCentroid
            x(f)=double(Spots.Fits(f).xDoG(Spots.Fits(f).z==Spots.Fits(f).brightestZ));
            y(f)=double(Spots.Fits(f).yDoG(Spots.Fits(f).z==Spots.Fits(f).brightestZ));
            z(f)=double(Spots.Fits(f).brightestZ);
        else
            x(f)=double(round(Spots.Fits(f).GaussPos(1)));
            y(f)=double(round(Spots.Fits(f).GaussPos(2)));
            z(f)=double(round(Spots.Fits(f).GaussPos(3)));
        end
    end
else
    x=[];
    y=[];
    z=[];
end

