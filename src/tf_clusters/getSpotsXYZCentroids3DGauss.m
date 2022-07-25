function [x, y, z, xyzGauss3D] = getSpotsXYZCentroids3DGauss(SpotsFrame)

% Return the X, Y, and Z positions, in pixels, of all the Spots in this 
% frame, using centroids from the 3D Gaussian fits

SpotsFits = SpotsFrame.Fits;

if isempty(SpotsFits)
    x = [];
    y = [];
    z = [];
    xyzGauss3D = [];

elseif ~isfield(SpotsFits,'GaussPos3D')
    error('No 3D Gaussian fit data detected. Have you run segementSpots with the ''fit3D'' or ''fit3Donly'' option?')

else
    nSpots = length(SpotsFits);
    
    x = nan(1, nSpots);
    y = nan(1, nSpots);
    z = nan(1, nSpots);
    xyzGauss3D = nan(nSpots, 3);
    
    for s = 1:nSpots
        gaussPos3D = SpotsFits(s).GaussPos3D;
        
        % Spots.mat stores all spot positions in units of pixels
        x(s) = gaussPos3D(1);
        y(s) = gaussPos3D(2);
        z(s) = gaussPos3D(3);
        xyzGauss3D(s,1) = gaussPos3D(1);
        xyzGauss3D(s,2) = gaussPos3D(2);
        xyzGauss3D(s,3) = gaussPos3D(3);
    end
end

