function [x,y,z,f]=SpotsXYZ(Spots, varargin)

useGauss3DCentroid = false;
if ~isempty(varargin)
    useGauss3DCentroid = varargin{1};
end

%Return the X and Y coordinate of the brightest Z of each spot in the
%Spots structure

nFits = length(Spots.Fits);

if ~isempty(Spots.Fits)
    
    x = nan(1, nFits);
    y = nan(1, nFits);
    z = nan(1, nFits);
    f = nan(1, nFits);
    for frame = 1:nFits
        
        spotsFrame = Spots.Fits(frame);
        
        if ~useGauss3DCentroid
            brightestZ = spotsFrame.brightestZ;
            brightestZIndex = spotsFrame.z == brightestZ;
            x(frame)=double(spotsFrame.xFit(brightestZIndex));
            y(frame)=double(spotsFrame.yFit(brightestZIndex));
            f(frame)=double(spotsFrame.FixedAreaIntensity(brightestZIndex));
            z(frame)=double(brightestZ);            
        else
            x(frame)=double(spotsFrame.GaussPos3D(2));
            y(frame)=double(spotsFrame.GaussPos3D(1));
            z(frame)=double(spotsFrame.GaussPos3D(3));
            f(frame)=double(spotsFrame.gauss3DIntensity);
        end
    end
else
    x=[];
    y=[];
    z=[];
    f=[];
end

