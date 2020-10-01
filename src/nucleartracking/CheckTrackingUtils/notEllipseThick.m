function h=notEllipse(ra,rb,ang,x0,y0,C,Nb, ax)
% Ellipse adds ellipses to the current plot
%
% MODIFIED BY AR 7/9/18: NOW REQUIRES A FULL 8 INPUTS FOR EVERY FUNCTION
% CALL AND THE AXES HANDLE MUST BE SPECIFIED AS THE LAST INPUT. THE REST OF THE
% DOCUMENTATION DOES NOT REFLECT THIS CHANGE. 
%
% ELLIPSE(ra,rb,ang,x0,y0) adds an ellipse with semimajor axis of ra,
% a semimajor axis of radius rb, a semimajor axis of ang, centered at
% the point x0,y0.
%
% The length of ra, rb, and ang should be the same. 
% If ra is a vector of length L and x0,y0 scalars, L ellipses
% are added at point x0,y0.
% If ra is a scalar and x0,y0 vectors of length M, M ellipse are with the same 
% radii are added at the points x0,y0.
% If ra, x0, y0 are vectors of the same length L=M, M ellipses are added.
% If ra is a vector of length L and x0, y0 are  vectors of length
% M~=L, L*M ellipses are added, at each point x0,y0, L ellipses of radius ra.
%
% ELLIPSE(ra,rb,ang,x0,y0,C)
% adds ellipses of color C. C may be a string ('r','b',...) or the RGB value. 
% If no color is specified, it makes automatic use of the colors specified by 
% the axes ColorOrder property. For several circles C may be a vector.
%
% ELLIPSE(ra,rb,ang,x0,y0,C,Nb), Nb specifies the number of points
% used to draw the ellipse. The default value is 300. Nb may be used
% for each ellipse individually.
%
% h=ELLIPSE(...) returns the handles to the ellipses.
%
% as a sample of how ellipse works, the following produces a red ellipse
% tipped up at a 45 deg axis from the x axis
% ellipse(1,2,pi/8,1,1,'r')
%
% note that if ra=rb, ELLIPSE plots a circle
%

% written by D.G. Long, Brigham Young University, based on the
% CIRCLES.m original 
% written by Peter Blattner, Institute of Microtechnology, University of 
% Neuchatel, Switzerland, blattner@imt.unine.ch

 
if nargin<6,
  C=[];
end

if nargin<7,
  Nb=[];
end

if isempty(Nb),Nb=300;end;
if isempty(C),C='b';end;


% work on the variable sizes

x0=x0(:);
y0=y0(:);
ra=ra(:);
rb=rb(:);
ang=ang(:);
Nb=Nb(:);

if isstr(C),C=C(:);end;

% how many inscribed elllipses are plotted

if length(ra)~=length(x0)
  maxk=length(ra)*length(x0);
else
  maxk=length(ra);
end;

% drawing loop
the = (0:Nb-1).*2*pi/(Nb-1);
the(1) = 0;
the(end) = 2*pi;
costhe = cos(the);
sinthe = sin(the);
CD = cos(ang);
SD = sin(ang);
o = ones(1,Nb);
Ax = [diag(CD) -diag(SD) x0];
Ay = [diag(SD) diag(CD) y0];
B = [ra*costhe;rb*sinthe; o];
X = Ax*B;
Y = Ay*B;

h =line(ax, X',Y','color',C, 'LineWidth', 2);

end
