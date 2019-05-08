function G = gauss3D(sigma, varargin)
 
    if isempty(varargin)
        kernelScale = 3;
        kernelBound = kernelScale*sigma + 1;
    else
        customSize = varargin{1};
        kernelBound = customSize;
    end
    x = -kernelBound:1:kernelBound;
    y = x;
    z = x;
    [x_, y_,z_] = meshgrid(x,y,z);
    A = (1/(sigma*sqrt(2*pi)))^3;
    G = A*exp((1/2)*(-(x_.^2)/(2*sigma^2) - (y_.^2)/(2*sigma^2) - (z_.^2)/(2*sigma^2)));
end
