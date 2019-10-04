function snip3D = simulate3DGaussGeneral(dim_vec, params)                                                       

% extract dim info
xDim = dim_vec(2);
yDim = dim_vec(1);
zDim = dim_vec(3);

% define grid ref arrays
[mesh_y,mesh_x, mesh_z] = meshgrid(1:yDim, 1:xDim, 1:zDim); 

% define helper arrays
mesh_x3 = @(params) reshape(mesh_x-params(3),1,1,[]);
mesh_y3 = @(params) reshape(mesh_y-params(2),1,1,[]);
mesh_z3 = @(params) reshape(mesh_z-params(4),1,1,[]);
xyz_vec = @(params) vertcat(mesh_x3(params),mesh_y3(params),mesh_z3(params));
xyz_array = @(params) repmat(xyz_vec(params),1,3,1);

% define covariance matrix
chol_mat = @(params) [params(5),    0,    0;
                      params(8),    params(6),    0;
                      params(9),    params(10),    params(7)];
               
% cov_mat = @(params) [params(5)^2,    params(8)*params(5)*params(6),    params(9)*params(5)*params(7);
%                    params(8)*params(5)*params(6),    params(6)^2,    params(10)*params(6)*params(7);
%                    params(9)*params(5)*params(7),    params(10)*params(6)*params(7),    params(7)^2];
cov_mat = @(params) chol_mat(params)*transpose(chol_mat(params));               
% define function to invert covariance array and expand to be size of mesh
cov_inv_array = @(params) repmat(inv(cov_mat(params)),1,1,numel(mesh_x));

% define Gauss function
single3DGaussian = @(params) params(1).*...
    exp(-.5*reshape(sum(reshape(sum(xyz_array(params).*cov_inv_array(params),1),3,1,[]).*...
    xyz_vec(params),1),size(mesh_x,1),size(mesh_x,2),size(mesh_x,3)));

snip3D = single3DGaussian(params);