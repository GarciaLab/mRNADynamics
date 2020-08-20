clear
close all

snip3D = rand(15,15,15);

% define grid ref arrays
[mesh_y,mesh_x, mesh_z] = meshgrid(1:size(snip3D,2), 1:size(snip3D,1), 1:size(snip3D, 3));

% initialize parameters
initial_parameters =[.1, 8,8, 8,1,1,1,0,0,0,.04];
test_parameters =[.1, 8,9, 7,2,1,2,.5,.3,-.4,.07];
ub_vec = [10,15,15,15,4,4,4,1,1,1,1];
lb_vec = [0,1,1,1,.5,.5,.5,-1,-1,-1,0];
% params = initial_parameters;    


mesh_x3 = @(params) reshape(mesh_x-params(3),1,1,[]);
mesh_y3 = @(params) reshape(mesh_y-params(2),1,1,[]);
mesh_z3 = @(params) reshape(mesh_z-params(4),1,1,[]);
xyz_vec = @(params) vertcat(mesh_x3(params),mesh_y3(params),mesh_z3(params));
xyz_array = @(params) repmat(xyz_vec(params),1,3,1);

% define covariance matrix
cov_mat = @(params) [params(5)^2 params(8)*params(5)*params(6) params(9)*params(5)*params(7);
                   params(8)*params(5)*params(6) params(6)^2 params(10)*params(6)*params(7);
                   params(9)*params(5)*params(7) params(10)*params(6)*params(7) params(7)^2];
               
cov_array = @(params) repmat(inv(cov_mat(params)),1,1,numel(mesh_x));

% define Gauss function
single3DGaussian = @(params) params(1).*...
exp(-.5*reshape(sum(reshape(sum(xyz_array(params).*cov_array(params),1),3,1,[]).*xyz_vec(params),1),size(mesh_x,1),size(mesh_x,2),size(mesh_x,3)));

snip3D = single3DGaussian(test_parameters) + initial_parameters(end);
% define objective function
single3DObjective = @(params) single3DGaussian(params) + params(11) - double(snip3D);     
% attempt to fit
meh = lsqnonlin(single3DObjective,initial_parameters,lb_vec,ub_vec);

