function [x,y,z]=getSpotsXYZ(Spots, varargin)

useGauss3DCentroid = false;
if ~isempty(varargin)
    useGauss3DCentroid = true;
end

%Return the X and Y coordinate of the brightest Z of each spot in the
%Spots structure

nFrames = length(Spots.Fits);

if ~isempty(Spots.Fits)
    
    x = nan(1, nFrames);
    y = nan(1, nFrames);
    z = nan(1, nFrames);
    
    for frame = 1:nFrames
        
        spotsFrame = Spots.Fits(frame);
        
        if ~useGauss3DCentroid
            brightestZ = spotsFrame.brightestZ;
            brightestZIndex = spotsFrame.z == brightestZ;
            x(frame)=double(spotsFrame.xDoG(brightestZIndex));
            y(frame)=double(spotsFrame.yDoG(brightestZIndex));
            z(frame)=double(brightestZ);
        else
            x(frame)=double(round(spotsFrame.GaussPos(1)));
            y(frame)=double(round(spotsFrame.GaussPos(2)));
            z(frame)=double(round(spotsFrame.GaussPos(3)));
        end
    end
else
    x=[];
    y=[];
    z=[];
end

