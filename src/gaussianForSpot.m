function singleGaussian = gaussianForSpot(snippet)

    %fits: [amplitude, x position, x width, y position, y width, offset, angle] 

    [y,x] = meshgrid(1:size(snippet,2), 1:size(snippet,1));

    singleGaussian = @(params) (params(1).*...
    exp(-(...
    (((cos(params(7)))^2 / (2*params(3)^2) ) + ((sin(params(7)))^2 / (2*params(5)^2)))  .* (x-params(2)).^2 ...
    - 2*((-sin(2*params(7)) / (4*params(3)^2) ) + (sin(2*params(7)) / (4*params(5)^2))) .* (x-params(2)).*(y-params(4))...
    + (((sin(params(7)))^2 / (2*params(3)^2) ) + ((cos(params(7)))^2 / (2*params(5)^2))).* (y-params(4)).^2 ...
        )))...
    + params(6) - double(snippet);

    %linear version
    
    %lnF = (lnA -a*xo^2 - c*yo^2 - 2*b*xo*yo) 
    % - a*x^2
    % - c*y^2
    % + (2*a*xo + 2*b*yo)x
    % + (2*b*xo + 2*c*yo)y
    % + (-2b)x*y
    %so a = -1*params(1)
    % c = -1*params(6)
    % b = -(1/2) * params(4)
    % xo = (b*params(5) - c*params(3))/ (2*b^2 - 2*a*c)
    % yo = (params(5)/(2*b)) - (a/b)*xo
    % A = exp (params(1) + a*xo^2+c*yo^2+2*b*xo*yo) 
    
    logGaussian = @(params) params(1) + params(2)*x^2 + params(3)*x+...
        params(4)*x*y + params(5)*y + params(6)*y^2;
     

     
end