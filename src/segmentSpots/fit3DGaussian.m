function [GaussFit, FitDeltas, GaussIntegral, GaussIntegralSE] = fit3DGaussian(snip3D,PixelSize,varargin)
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
    % take a guess at sigma
    sigma_guess = 200/PixelSize;
    % define initial parameters
    xDim = size(snip3D,1);
    yDim = size(snip3D,2);
    zDim = size(snip3D,3);
    % initialize parameters
    initial_parameters =[max(snip3D(:)), ceil(yDim/2),ceil(xDim/2), ceil(zDim/2),sigma_guess,sigma_guess,sigma_guess,.1,.1,.1,prctile(snip3D(:),25)];
    
    % initialize upper and lower parameter bounds
    ub_vec = [Inf,yDim,xDim,zDim,xDim/8,xDim/8,xDim/8,xDim/8,xDim/8,xDim/8,Inf];
    lb_vec = [0,1,1,1,0,0,0,0,0,0,0];
    
    % check for additional arguments
    for i = 1:(numel(varargin)-1)  
        if ischar(varargin{i})
            eval([varargin{i} '= varargin{i+1};']);        
        end
    end
       
    % define objective function
    dim_vec = [yDim, xDim, zDim];    
    single3DObjective = @(params) simulate3DGaussGeneral(dim_vec, params) + params(end) - double(snip3D);     
    % attempt to fit
    options.Display = 'off';
    [GaussFit, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(single3DObjective,initial_parameters,lb_vec,ub_vec,options);
    
    % estimate confidence intervals for gauss integral numerically (faster
    % than symbolics)
    FitCI = nlparci(GaussFit,residual,'jacobian',jacobian);
    FitCI(FitCI<0) = 0;
    FitDeltas = diff(FitCI')' / 2 / 1.96;
    paramValues = normrnd(0,1,size(FitDeltas,1),100).*FitDeltas + GaussFit';
    paramValues(paramValues<0) = realmin;

    chol_mat = @(params) [params(5),    0,    0;
                          params(8),    params(6),    0;
                          params(9),    params(10),    params(7)];
    cov_mat = @(params) chol_mat(params)*transpose(chol_mat(params));   
    gauss_int = @(params) params(1)*(2*pi)^1.5 *sqrt(det(cov_mat(params)));
    gauss_int_vec = NaN(1,100);
    for i = 1:size(paramValues,2)
        gauss_int_vec(i) = gauss_int(paramValues(:,i));
    end
    % extract values to report
    GaussIntegral = gauss_int(GaussFit);
    GaussIntegralSE = nanstd(gauss_int_vec);
    
    % perform some QC checks on fit
    
%     x_sigma = 2;
%     y_sigma = 2;
%     z_sigma = 2;
%     % estimate raw fluorescence
%     [x_ref_vol, y_ref_vol, z_ref_vol] = meshgrid(1:size(snip3D,1),1:size(snip3D,2),1:size(snip3D,3));
%     dx = x_ref_vol - GaussFit(2);
%     dy = y_ref_vol - GaussFit(3);
%     dz = z_ref_vol - GaussFit(4);
%     
%     wt_array = exp(-.5*((dx/x_sigma).^2 + (dy/y_sigma).^2 + (dz/z_sigma)^2));
%     wt_array = wt_array / sum(wt_array(:));
%     
%     GaussIntegralRaw = sum(wt_array(:).*snip3D(:)) - sum(numel(snip3D)*GaussFit(end)*(wt_array(:)));
end