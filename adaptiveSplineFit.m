function optFit = adaptiveSplineFit(x,y,n, R)
% find the best spline fit with n breaks
if nargin<4
    R=1e-5; % "robust" mode off
elseif strcmpi(R,'r')
    R=0.5;
end

% position of the breaks: at the endpoints of the range, and n-2 somewhere
% in between
a = min(x);
b = max(x);
breaks = linspace(a, b, n);
x0 = breaks(2:end-1);

lb = a*ones(size(x0));
ub = b*ones(size(x0));

options = optimset('Display','off');
optParams = lsqnonlin(@(breaks)y-spfit(x, y, a, b, breaks), ...
    x0, lb, ub, options);

% force zero-derivative boundary conditions
xc = [a b];
yc = [0 0];
cc = [0 0; 1 1];
con = struct('xc',xc,'yc',yc,'cc',cc);
optFit = splinefit(x,y,[a sort(optParams) b], 4, R, con);
end

function yfit = spfit(x, y, a, b, params)
    breaks = [a params b];
    pp = splinefit(x,y,sort(breaks), 4);
    yfit = ppval(pp, x);
    % force the break points to be at least 0.01 apart
    penalty = -sum(min(diff(breaks)-0.01, 0)); % if all steps are larger than 0.01, this is zero
    yfit = yfit + 100*penalty;
end