function snip3D = simulate3DGaussianRhoSimple(dimensionVector, params, inverseCovarianceMatrix)                                                       

% extract dim info
xDim = dimensionVector(2);
yDim = dimensionVector(1);
zDim = dimensionVector(3);

% define grid ref arrays
[mesh_y,mesh_x, mesh_z] = meshgrid(1:yDim, 1:xDim, 1:zDim); 

% define helper arrays
mesh_x3 = @(params) reshape(mesh_x-params(3),1,1,[]);
mesh_y3 = @(params) reshape(mesh_y-params(2),1,1,[]);
mesh_z3 = @(params) reshape(mesh_z-params(4),1,1,[]);
xyzVector = @(params) vertcat(mesh_x3(params),mesh_y3(params),mesh_z3(params));
xyzArray = @(params) repmat(xyzVector(params),1,3,1);

% define covariance matrix               
% cov_mat = @(params) [params(5)^2,    params(8)*params(5)*params(6),    params(9)*params(5)*params(7);
%                    params(8)*params(5)*params(6),    params(6)^2,    params(10)*params(6)*params(7);
%                    params(9)*params(5)*params(7),    params(10)*params(6)*params(7),    params(7)^2];
               
% define function to invert covariance array and expand to be size of mesh
covarianceInverseArray = repmat(inverseCovarianceMatrix,1,1,numel(mesh_x));

% define Gauss function
single3DGaussian = @(params) params(1).*...
    exp(-.5*reshape(sum(reshape(sum(xyzArray(params).*covarianceInverseArray,1),3,1,[]).*...
    xyzVector(params),1),size(mesh_x,1),size(mesh_x,2),size(mesh_x,3)));

snip3D = single3DGaussian(params);