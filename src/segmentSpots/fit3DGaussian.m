function [GaussFit, SpotIntegral, GaussCI, IntegralCI, OffsetCI] = fit3DGaussian(snip3D,PixelSize,varargin)
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
    initial_parameters =[max(snip3D(:)), ceil(yDim/2),ceil(xDim/2), ceil(zDim/2),sigma_guess,sigma_guess,prctile(snip3D(:),25)];
    
    % initialize upper and lower parameter bounds
    ub_vec = [Inf,yDim,xDim,zDim,xDim,zDim,Inf];
    lb_vec = [0,1,1,1,.5,.5,0];
    
    % check for additional arguments
    for i = 1:(numel(varargin)-1)  
        if ischar(varargin{i})
            eval([varargin{i} '= varargin{i+1};']);        
        end
    end
       
    % define objective function
    dim_vec = [yDim, xDim, zDim];    
    single3DObjective = @(params) simulate3DGaussSymmetric(dim_vec, params) - double(snip3D);     
    % attempt to fit
    [GaussFit, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(single3DObjective,initial_parameters,lb_vec,ub_vec);
    FitCI = nlparci(GaussFit,residual,'jacobian',jacobian);
    GaussCI = FitCI(1:6,:);
    Deltas = diff(gaussCI');
    OffsetCI = FitCI(end,:);
    % perform error propagation to obtain error bounds on integral   
    sigmaXY = GaussFit(5);
    sigmaZ = GaussFit(6);
    amplitude = GaussFit(1);
    
    SpotIntegral = amplitude * (2*pi)^1.5 * sigmaXY^2 * sigmaZ;   
    IntError = (2*pi)^1.5 * sqrt((Deltas(1)*sigmaXY^2 * sigmaZ)^2 + (Deltas(6)*sigmaXY^2 * amplitude)^2 ...
        + (2*Deltas(5)*amplitude * sigmaXY^2 * sigmaZ)^2);    
    IntegralCI =  [SpotIntegral - IntError/2, SpotIntegral + IntError/2];
end