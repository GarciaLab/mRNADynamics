%   Function: compute_directellipse_estimates
%
%   This function is a wrapper for the numerically stable direct ellipse
%   fit due to 
%
%   R. Halif and J. Flusser
%   "Numerically stable direct least squares fitting of ellipses"
%   Proc. 6th International Conference in Central Europe on Computer 
%   Graphics and Visualization. WSCG '98 Czech Republic,125--132, feb, 1998
%
%   which is a modificaiton of
%
%   A. W. Fitzgibbon, M. Pilu, R. B. Fisher
%   "Direct Least Squares Fitting of Ellipses"
%   IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
%
%   It first shifts all the data points to a new
%   coordinate system so that the origin of the coordinate system is at the
%   center of the data points, and then scales all data points so that they
%   lie more or less within a unit box. It is within this transformed
%   coordinate system that the ellipse is estimated. The resulting ellipse
%   parameters are then transformed back to the original data space. 
%
%   Parameters:
%
%      dataPts    - a Nx2 matrix where N is the number of data points	
%
%   Returns: 
%
%     a length-6 vector [a b c d e f] representing the parameters of the
%     equation
%    
%     a x^2 + b x y + c y^2 + d x + e y + f = 0
%
%     with the additional result that b^2 - 4 a c < 0.
%
%   See Also: 
%
%    fast_guaranteed_ellipse_estimate
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 27/3/2014 
function [theta]  = compute_directellipse_estimates(dataPts)


 nPts = length(dataPts);

% scale and translate data points so that they lie inside a unit box
[normalizedPoints, T] = normalize_data_isotropically(dataPts);
normalizedPoints = [ normalizedPoints, ones( nPts,1 ) ];

theta = direct_ellipse_fit(normalizedPoints'); 
theta = theta / norm(theta);

a = theta(1); b = theta(2) ; c = theta(3) ;
d = theta(4) ; e = theta(5); f = theta(6);
C = [a b/2 d/2 ; b/2 c e/2; d/2 e/2 f];

% denormalise C
C = T'*C*T;
aa = C(1,1);
bb = C(1,2)*2;
dd = C(1,3)*2;
cc = C(2,2);
ee = C(2,3)*2;
ff = C(3,3);
theta = [aa bb cc dd ee ff]';
theta = theta / norm(theta);

end

