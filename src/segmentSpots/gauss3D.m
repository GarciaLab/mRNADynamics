function G = gauss3D(sigmaXY, varargin)

kernelScale = 3;
kernelBound = kernelScale*sigmaXY + 1;
sigmaZ = sigmaXY;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'kernelSize')
        kernelBound = varargin{i+1};
    elseif strcmpi(varargin{i}, 'sigmaZ')
        sigmaZ = varargin{i+1};
    end
end

x = -kernelBound:1:kernelBound;
y = x;
z = x;
[x_, y_,z_] = meshgrid(x,y,z);
A = (1/(sigmaXY*sqrt(2*pi)))^3;
G = A*exp((1/2)*(-(x_.^2)/(2*sigmaXY^2) - (y_.^2)/(2*sigmaXY^2) - (z_.^2)/(2*sigmaZ^2)));

end
