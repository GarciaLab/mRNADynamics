function G = gauss3D(sigma)

    C = 3;
    M = C*sigma + 1;
    x = -M:1:M;
    y = x;
    z = x;
    sigma = 3;
    [x_, y_,z_] = meshgrid(x,y,z);
    A = 1/sqrt(2*pi*sigma^2)^3;
    G = A*exp(-x_.^2/(2*sigma^2) + -y_.^2/(2*sigma^2) + -z_.^2/(2*sigma^2));
    
end
