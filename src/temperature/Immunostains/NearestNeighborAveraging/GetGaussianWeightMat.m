function GaussianWeights = GetGaussianWeightMat(x, y, sigma, width)
%%
if ~exist('width', 'var')
    width = [];
end
Nx = length(x);
Ny = length(y);

GaussianWeights = zeros(Nx,Ny);
for i = 1:Nx
    for j =1:Ny
        GaussianWeights(i,j) = exp(-((x(i)-y(j))^2)/(2*sigma^2));
    end
end
if ~isempty(width)
    MinWeight = exp(-(width^2)/(2*sigma^2));
    GaussianWeights(GaussianWeights<MinWeight) = 0;
end

