function [x,y,z] = generateRandomSpots(nSpots,ellipseCentX,ellipseCentY,ellipseCentZ,nuclearRadius)
% Generates N random points in a circle.
% RAND_CIRC(N) generates N random points in the unit circle at (0,0).
% RAND_CIRC(N,x,y,r) generates N random points in a circle with radius r 
% and center at (x,y).
if nargin<2
   ellipseCentX = 0;
   ellipseCentY = 0;
   ellipseCentZ = 0;
   nuclearRadius = 1;
end

rvals = 2*rand(nSpots,1)-1;
elevation = asin(rvals);
azimuth = 2*pi*rand(nSpots,1);
radii = nuclearRadius*(rand(nSpots,1).^(1/3));

[xCoord,yCoord,zCoord] = sph2cart(azimuth,elevation,radii);

x = xCoord + ellipseCentX;
y = yCoord + ellipseCentY;
z = zCoord + ellipseCentZ;