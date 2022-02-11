function snip3D = simulate3DGaussSymmetric(dimensionVector, params)

% extract dim info
xDim = dimensionVector(2);
yDim = dimensionVector(1);
zDim = dimensionVector(3);

% define grid ref arrays
[mesh_y, mesh_x, mesh_z] = meshgrid(1:yDim, 1:xDim, 1:zDim); 

% define helper arrays
mesh_x3 = @(params) reshape(mesh_x-params(5),1,1,[]);
mesh_y3 = @(params) reshape(mesh_y-params(4),1,1,[]);
mesh_z3 = @(params) reshape(mesh_z-params(6),1,1,[]);
xyzVector = @(params) vertcat(mesh_x3(params),mesh_y3(params),mesh_z3(params));
xyzArray = @(params) repmat(xyzVector(params),1,3,1);

% define covariance matrix
covarianceMatrix = @(params) [params(2)^2       0              0;
                                  0      params(2)^2           0;
                                  0            0        params(3)^2];
               
% define function to invert covariance array and expand to be size of mesh
covarianceArray = @(params) repmat(inv(covarianceMatrix(params)),1,1,numel(mesh_x));

% define Gauss function
single3DGaussian = @(params) params(1).*...
    exp(-.5*reshape(sum(reshape(sum(xyzArray(params).*covarianceArray(params),1),3,1,[]).*...
    xyzVector(params),1),size(mesh_x,1),size(mesh_x,2),size(mesh_x,3)));

snip3D = single3DGaussian(params);