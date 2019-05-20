function singleGaussian = gaussianForSpot(y, x, snippet)

    %fits: [amplitude, x position, x width, y position, y width, offset, angle, linear offset x, linear offset y] 

    singleGaussian = @(params) (params(1).*...
    exp(-(...
    (((cos(params(7)))^2 / (2*params(3)^2) ) + ((sin(params(7)))^2 / (2*params(5)^2)))  .* (x-params(2)).^2 ...
    - 2*((-sin(2*params(7)) / (4*params(3)^2) ) + (sin(2*params(7)) / (4*params(5)^2))) .* (x-params(2)).*(y-params(4))...
    + (((sin(params(7)))^2 / (2*params(3)^2) ) + ((cos(params(7)))^2 / (2*params(5)^2))).* (y-params(4)).^2 ...
        )))...
    + params(6) + params(8).*x + params(9).*y - double(snippet);
     
end