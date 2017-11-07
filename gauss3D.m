function G = gauss3D(sigma, customSize)

    if isempty(customSize)
        C = 3;
        M = C*sigma + 1;
    else
        M = customSize{1};
    end
    x = -M:1:M;
    y = x;
    z = x;
    [x_, y_,z_] = meshgrid(x,y,z);
    A = (1/(sigma*sqrt(2*pi)))^3;
    G = A*exp((1/2)*(-(x_.^2)/(2*sigma^2) - (y_.^2)/(2*sigma^2) - (z_.^2)/(2*sigma^2)));
end
