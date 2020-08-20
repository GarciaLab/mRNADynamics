% INPUT ARGUMENTS:
% snip3D: 3D array containing spot to fit. Should contain only one spot

% OPTIONS:
% initial_parameters: 1 x 11 vector containing values from which to
%                     initialize inferemce
% ub_vec: 1 x 11 vector specifying upper bounds for inference parameters
% lb_vec: 1 x 11 vector specifying lower bounds

% RETURNS:
% GaussFit: 1 x 11 vector of inferred parameter values
% Parameter identity is as follows: 
%        (1) amplitude of gaussian
%        (2-4) y,x,and z center positions
%        (5-7) y,x,and z sigma values
%        (8-10) y,x,and z sigma covariance coefficients
%        (11) inferred background offset

function GaussFit = fit3DGaussian(snip3D,varargin)

% define initial parameters
xDim = size(snip3D,1);
yDim = size(snip3D,2);
zDim = size(snip3D,3);
% initialize parameters
initial_parameters =[.1, ceil(yDim/2),ceil(xDim/2), ceil(zDim/2),1,1,1,0,0,0,.04];

% initialize upper and lower parameter bounds
ub_vec = [max(snip3D(:)),yDim,xDim,zDim,yDim,xDim,zDim,1,1,1,median(snip3D(:))];
lb_vec = [0,1,1,1,.5,.5,.5,-1,-1,-1,0];

% check for additional arguments
for i = 1:(numel(varargin)-1)  
    if ischar(varargin{i})
        eval([varargin{i} '= varargin{i+1};']);        
    end
end

% define grid ref arrays
[mesh_y,mesh_x, mesh_z] = meshgrid(1:size(snip3D,2), 1:size(snip3D,1), 1:size(snip3D, 3));  

% define helper arrays
mesh_x3 = @(params) reshape(mesh_x-params(3),1,1,[]);
mesh_y3 = @(params) reshape(mesh_y-params(2),1,1,[]);
mesh_z3 = @(params) reshape(mesh_z-params(4),1,1,[]);
xyz_vec = @(params) vertcat(mesh_x3(params),mesh_y3(params),mesh_z3(params));
xyz_array = @(params) repmat(xyz_vec(params),1,3,1);

% define covariance matrix
cov_mat = @(params) [params(5)^2 params(8)*params(5)*params(6) params(9)*params(5)*params(7);
                   params(8)*params(5)*params(6) params(6)^2 params(10)*params(6)*params(7);
                   params(9)*params(5)*params(7) params(10)*params(6)*params(7) params(7)^2];
 
% define function to invert covariance array and expand to be size of mesh
cov_array = @(params) repmat(inv(cov_mat(params)),1,1,numel(mesh_x));

% define Gauss function
single3DGaussian = @(params) params(1).*...
    exp(-.5*reshape(sum(reshape(sum(xyz_array(params).*cov_array(params),1),3,1,[]).*...
    xyz_vec(params),1),size(mesh_x,1),size(mesh_x,2),size(mesh_x,3)));

% snip3D = single3DGaussian(test_parameters) + initial_parameters(end) + rand(size(snip3D))/50;
% define objective function
single3DObjective = @(params) single3DGaussian(params) + params(11) - double(snip3D);     
% attempt to fit
GaussFit = lsqnonlin(single3DObjective,initial_parameters,lb_vec,ub_vec);