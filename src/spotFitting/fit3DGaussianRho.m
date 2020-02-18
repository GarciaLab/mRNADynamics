function [GaussFit, FitDeltas, GaussIntegral, GaussIntegralSE, GaussIntegralRaw, np_flag] = ...
                                                    fit3DGaussianRho(snip3D,PixelDims,varargin)
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
    % generate ref arrays
    [x_ref_vol, y_ref_vol, z_ref_vol] = meshgrid(1:size(snip3D,1),1:size(snip3D,2),1:size(snip3D,3));
    % take a guess at sigma
    sigma_guess = 200/PixelDims(1);
    % set dimensions of integration volume for raw estimate
    xy_sigma = 300 / PixelDims(1);
    z_sigma = 750 / PixelDims(2);
    % define initial parameters
    xDim = size(snip3D,1);
    yDim = size(snip3D,2);
    zDim = size(snip3D,3);
    intMax = max(snip3D(:));
    int25 = prctile(snip3D(:),25);
    % initialize parameters
    initial_parameters =[intMax, ceil(yDim/2),ceil(xDim/2), ceil(zDim/2),sigma_guess,sigma_guess,sigma_guess,0,0,0,int25,0,0,0];
    
    % initialize upper and lower parameter bounds
    ub_vec = [intMax*1.5,yDim,xDim,zDim,xDim/8,xDim/8,xDim/8,1,1,1,intMax,intMax/yDim,intMax/xDim,intMax/zDim];
    lb_vec = [0,1,1,1,1,1,1,-1,-1,-1,0,-intMax/yDim,-intMax/xDim,-intMax/zDim];
    
    % check for additional arguments
    for i = 1:(numel(varargin)-1)  
        if ischar(varargin{i})
            eval([varargin{i} '= varargin{i+1};']);        
        end
    end
       
    % define objective function
    dim_vec = [yDim, xDim, zDim];    
    single3DObjective = @(params) simulate3DGaussianRho(dim_vec, params) + params(11) + ...
            params(12)*x_ref_vol + params(13)*y_ref_vol + params(14)*z_ref_vol - double(snip3D);     
    % attempt to fit
    options.Display = 'off';
    [GaussFit, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(single3DObjective,initial_parameters,lb_vec,ub_vec,options);
    % check that inferred covariance matrix is positive semi-definite
    cov_mat = @(params) [params(5)^2,    params(8)*params(5)*params(6),    params(9)*params(5)*params(7);
                   params(8)*params(5)*params(6),    params(6)^2,    params(10)*params(6)*params(7);
                   params(9)*params(5)*params(7),    params(10)*params(6)*params(7),    params(7)^2];
               
    [~, np_flag] = chol(cov_mat(GaussFit));
    if np_flag 
        % find nearest SPD matrix
        cov_new = nearestSPD(cov_mat(GaussFit));
        % update parameters
        GaussFitNew = NaN(size(GaussFit));
        
        GaussFitNew(5:7) = sqrt(diag(cov_new));
        off_diag_vec = cov_new([2 3 6]);
        GaussFitNew(8:10) = off_diag_vec ./ GaussFitNew([5 5 6]) ./ GaussFitNew([6 7 7]);
        % conduct secondary inference to adjust other params
        single3DObjectiveSimp = @(params) simulate3DGaussianRhoSimple(dim_vec, params,inv(cov_new)) + params(11) + ...
            params(12)*x_ref_vol + params(13)*y_ref_vol + params(14)*z_ref_vol - double(snip3D);  
        [GaussFitNew([1:4 11:14]), ~, ~, ~, ~, ~, ~] = lsqnonlin(single3DObjectiveSimp,initial_parameters([1:4 11:14]),lb_vec([1:4 11:14]),ub_vec([1:4 11:14]),options);
        GaussFit = GaussFitNew;
    end
        
        
    % estimate confidence intervals for gauss integral numerically (faster
    % than symbolics)
    gauss_int = @(params) params(1)*(2*pi)^1.5 *sqrt(det(cov_mat(params)));
    GaussIntegralSE = NaN;
    GaussIntegral = gauss_int(GaussFit);  
    FitDeltas = NaN(size(GaussFit));
    if ~np_flag
        FitCI = nlparci(GaussFit,residual,'jacobian',jacobian);
        FitCI(FitCI<0) = 0;
        FitDeltas = diff(FitCI')' / 2 / 1.96;
        paramValues = normrnd(0,1,size(FitDeltas,1),100).*FitDeltas + GaussFit';
        paramValues(paramValues<0) = realmin;
        % estimate gauss integral
        gauss_int_vec = NaN(1,100);
        for i = 1:size(paramValues,2)
            gauss_int_vec(i) = gauss_int(paramValues(:,i));
        end
        GaussIntegralSE = nanstd(gauss_int_vec);
%         GaussIntegral = nanmean(gauss_int_vec);
    end                           
    bkg_inf = GaussFit(11) + GaussFit(12)*x_ref_vol + GaussFit(13)*y_ref_vol + GaussFit(14)*z_ref_vol ;
    % estimate raw fluorescence    
    dx = abs(x_ref_vol - GaussFit(2));
    dy = abs(y_ref_vol - GaussFit(3));
    dz = abs(z_ref_vol - GaussFit(4));
    wt_array = double((dx <= 1.96*xy_sigma).*(dy<=1.96*xy_sigma).*(dz<=1.96*z_sigma)); 
    GaussIntegralRaw = sum(wt_array(:).*double(snip3D(:))) - sum(bkg_inf(:).*wt_array(:));
end