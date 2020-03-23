function [GaussFit, FitDeltas, GaussIntegral, GaussIntegralSE, GaussIntegralRaw, np_flag, integralDimensions] = ...
                                                    fit3DGaussianRho(snip3D,PixelDims,varargin)
    % INPUT ARGUMENTS:
    % snip3D: 3D array containing spot to fit. Should contain only one spot
    
    % OPTIONS:
    % initial_parameters: 1 x 11 vector containing values from which to
    %                     initialize inferemce
    % upperBoundVector: 1 x 11 vector specifying upper bounds for inference parameters
    % lowerBoundVector: 1 x 11 vector specifying lower bounds
    
    % RETURNS:
    % GaussFit: 1 x 11 vector of inferred parameter values
    % Parameter identity is as follows: 
    %        (1) amplitude of gaussian
    %        (2-4) y,x,and z center positions
    %        (5-7) y,x,and z sigma values
    %        (8-10) y,x,and z sigma covariance coefficients
    %        (11) inferred background offset
    % generate grids
    [xMesh, yMesh, zMesh] = meshgrid(1:size(snip3D,1),1:size(snip3D,2),1:size(snip3D,3));
    % take a guess at sigma
    sigma_guess = 200/PixelDims(1);
    % set dimensions of integration volume for raw estimate
    xy_sigma = 230 / PixelDims(1);
    z_sigma = 620 / PixelDims(2);
    integralDimensions = 1.96*[xy_sigma z_sigma];
    % define initial parameters
    xDim = size(snip3D,1);
    yDim = size(snip3D,2);
    zDim = size(snip3D,3);
    intensityMax = max(snip3D(:));
    intensity25 = prctile(snip3D(:),25);
    % initialize parameters
    initial_parameters =[intensityMax, ceil(yDim/2),ceil(xDim/2), ceil(zDim/2),...
        sigma_guess,sigma_guess,sigma_guess,0,0,0,intensity25,0,0,0];
    
    % initialize upper and lower parameter bounds
    upperBoundVector = [intensityMax*1.5,yDim,xDim,zDim,xDim/8,...
        xDim/8,xDim/8,1,1,1,intensityMax,intensityMax/yDim,intensityMax/xDim,intensityMax/zDim];
    lowerBoundVector = [0,1,1,1,1,1,1,-1,-1,-1,0,-intensityMax/yDim,-intensityMax/xDim,-intensityMax/zDim];
    
    % check for additional arguments
    for i = 1:(numel(varargin)-1)  
        if ischar(varargin{i})
            eval([varargin{i} '= varargin{i+1};']);        
        end
    end
       
    % define objective function
    dimensionVector = [yDim, xDim, zDim];    
    single3DObjective = @(params) simulate3DGaussianRho(dimensionVector, params) + params(11) + ...
            params(12)*xMesh + params(13)*yMesh + params(14)*zMesh - double(snip3D);     
    % attempt to fit
    options.Display = 'off';
    [GaussFit, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(single3DObjective,double(initial_parameters),...
        double(lowerBoundVector),double(upperBoundVector),options);
    % check that inferred covariance matrix is positive semi-definite
    covarianceMatrix = @(params) [params(5)^2,    params(8)*params(5)*params(6),    params(9)*params(5)*params(7);
                   params(8)*params(5)*params(6),    params(6)^2,    params(10)*params(6)*params(7);
                   params(9)*params(5)*params(7),    params(10)*params(6)*params(7),    params(7)^2];
               
    [~, np_flag] = chol(covarianceMatrix(GaussFit));
    if np_flag 
        % find nearest SPD matrix
        covarianceMatrixNew = nearestSPD(covarianceMatrix(GaussFit));
        % update parameters
        GaussFitNew = NaN(size(GaussFit));
        
        GaussFitNew(5:7) = sqrt(diag(covarianceMatrixNew));
        offDiagonalVector = covarianceMatrixNew([2 3 6]);
        GaussFitNew(8:10) = offDiagonalVector ./ GaussFitNew([5 5 6]) ./ GaussFitNew([6 7 7]);
        % conduct secondary inference to adjust other params
        single3DObjectiveSimple = @(params) simulate3DGaussianRhoSimple(dimensionVector, params,inv(covarianceMatrixNew)) + params(11) + ...
            params(12)*xMesh + params(13)*yMesh + params(14)*zMesh - double(snip3D);  
        [GaussFitNew([1:4 11:14]), ~, ~, ~, ~, ~, ~] = lsqnonlin(single3DObjectiveSimple,...
            initial_parameters([1:4 11:14]),lowerBoundVector([1:4 11:14]),upperBoundVector([1:4 11:14]),options);
        GaussFit = GaussFitNew;
    end
        
        
    % estimate confidence intervals for gauss integral numerically (faster
    % than symbolics)
    gaussIntegralFunction = @(params) params(1)*(2*pi)^1.5 *sqrt(det(covarianceMatrix(params)));
    GaussIntegralSE = NaN;
    GaussIntegral = gaussIntegralFunction(GaussFit);  
    FitDeltas = NaN(size(GaussFit));
    if ~np_flag
        FitCI = nlparci(GaussFit,residual,'jacobian',jacobian);
        FitCI(FitCI<0) = 0;
        FitDeltas = diff(FitCI')' / 2 / 1.96;
        paramValues = normrnd(0,1,size(FitDeltas,1),100).*FitDeltas + GaussFit';
        paramValues(paramValues<0) = realmin;
        % estimate gauss integral
        gaussIntegralVector = NaN(1,100);
        for i = 1:size(paramValues,2)
            gaussIntegralVector(i) = gaussIntegralFunction(paramValues(:,i));
        end
        GaussIntegralSE = nanstd(gaussIntegralVector);
%         GaussIntegral = nanmean(gauss_int_vec);
    end
    
    backgroundInference = GaussFit(11) + GaussFit(12)*xMesh + GaussFit(13)*yMesh + GaussFit(14)*zMesh ;
    
    % estimate raw fluorescence    
    dx = abs(xMesh - GaussFit(2));
    dy = abs(yMesh - GaussFit(3));
    dz = abs(zMesh - GaussFit(4));
    weightArray = double((dx <= integralDimensions(1)).*(dy<=integralDimensions(1)).*(dz<=integralDimensions(2))); 
    GaussIntegralRaw = (sum(weightArray(:).*double(snip3D(:))) - sum(backgroundInference(:).*weightArray(:)));
end