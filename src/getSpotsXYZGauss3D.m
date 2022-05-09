function [x,y,z]=getSpotsXYZGauss3D(Spots)

%Return the X and Y coordinate of the brightest Z of each spot in the
%Spots structure
x=NaN(1,numel(Spots.Fits));
y=NaN(1,numel(Spots.Fits));
z=NaN(1,numel(Spots.Fits));
fits3D = 
if isfield(Spots.Fits,'gauss3DIntensity')        
    for i=1:length(Spots.Fits)
        GaussPos = Spots.Fits(i).GaussPos;
        fits3D = Spots.Fits(i).fits3D;        
        x(i)=GaussPos(1);
        y(i)=GaussPos(2);
        z(i)=GaussPos(3);        
    end
end

