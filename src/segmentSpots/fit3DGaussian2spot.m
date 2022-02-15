function [GaussParams1, GaussParams2, offset, GaussIntegralVector, centroid_mean, GaussSE1, GaussSE2, offsetSE,...
    GaussIntegralSEVector, centroid_se]...
        = fit3DGaussian2spot(snip3D,PixelSize, varargin)
    % INPUT ARGUMENTS:
    % snip3D: 3D array containing spot to fit. Should contain only one spot
    % PixelSize: size of pixels (in nm)
    % OPTIONS:
    % initial_parameters: 1 x 11 vector containing values from which to
    %                     initialize inferemce
    % ub_vec: 1 x 11 vector specifying upper bounds for inference parameters
    % lb_vec: 1 x 11 vector specifying lower bounds
    
    % RETURNS:
    % GaussFit: 1 x 11 vector of inferred parameter values
    % Parameter identity is as follows: 
    %        (1) amplitude of gaussian (spot 1)    
    %        (2-4) y,x,and z center positions (spot 1)    
    %        (5-7) y,x,and z sigma values (spot 1)  
    %        (8) amplitude of gaussian (spot 1)
    %        (9-11) dy, dx, and dz distances from spot 1 to spot 2
    %        (12-14) y,x,and z sigma values (spot 2)  
    %        (15) inferred background offset    
    
    sigma_guess = 200/PixelSize;
    % define initial parameters
    xDim = size(snip3D,1);
    yDim = size(snip3D,2);
    zDim = size(snip3D,3);
    % initialize parameters
    initial_parameters =[max(snip3D(:)), floor(yDim/2)-1,floor(xDim/2)-1, floor(zDim/2)-1,sigma_guess,sigma_guess...
        ,.1,ceil(yDim/2)+1,ceil(xDim/2)+1, ceil(zDim/2)+1,1,1,prctile(snip3D(:),10)];
    
    % initialize upper and lower parameter bounds
    upperBoundVector = [3*max(snip3D(:)),yDim,xDim,zDim,3*sigma_guess,3*sigma_guess,...
              3*max(snip3D(:)),yDim,xDim,zDim,3*sigma_guess,3*sigma_guess,...
              Inf];
    lowerBoundVector = [0,1,1,1,.5,.5,0,1,1,1,.1,.1,0];
    
    % check for additional arguments
    for i = 1:(numel(varargin)-1)  
        if ischar(varargin{i})
            eval([varargin{i} '= varargin{i+1};']);        
        end
    end
       
    % define objective function
    dimensionVector = [yDim, xDim, zDim];        
    double3DObjective = @(params) simulate3DGaussSymmetric(dimensionVector, params(1:6))...
        + simulate3DGaussSymmetric(dimensionVector, params(7:12)) ...
        + params(end) - double(snip3D);     
    % attempt to fit
    options.Display = 'off';
    [GaussFit, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(double3DObjective,...
        initial_parameters,lowerBoundVector,upperBoundVector,options);
    % store parameters
    GaussParams1 = GaussFit(1:6);
    GaussParams2 = GaussFit(7:12);        
    offset = GaussFit(end);
    
    % estimate error in integral calculations numeriucally
    FitCI = nlparci(GaussFit,residual,'jacobian',jacobian);
    FitCI(FitCI<0) = 0;
    FitDeltas = diff(FitCI')' / 2 / 1.96;
    paramValues = normrnd(0,1,size(FitDeltas,1),100).*FitDeltas + GaussFit';
    paramValues(paramValues<0) = realmin;
    
    covarianceMatrix = @(params) [params(5),    0,    0;
                          0,    params(5),    0;
                          0,    0,     params(6)];    
    gaussIntegral = @(params) params(1)*(2*pi)^1.5 *sqrt(det(covarianceMatrix(params)));    
    gaussIntegralVector1 = NaN(1,100);
    gaussIntegralVector = NaN(1,100);    
    centroidPositionVector = NaN(100,3); 
    for i = 1:size(paramValues,2)
        gaussIntegralVector1(i) = gaussIntegral(paramValues(1:6,i));
        gaussIntegralVector(i) = gaussIntegral(paramValues(7:12,i));        
        centroidPositionVector(i,:) = gaussIntegralVector1(i)*paramValues(2:4,i) + ...
            gaussIntegralVector1(i)*paramValues(8:10,i);
    end
    % extract values to report
    GaussIntegral1 = nanmean(gaussIntegralVector1);
    GaussIntegralSE1 = nanstd(gaussIntegralVector1);
    GaussIntegral2 = nanmean(gaussIntegralVector);
    GaussIntegralSE2 = nanstd(gaussIntegralVector);
    GaussIntegralTot = GaussIntegral1 + GaussIntegral2;
    GaussIntegralSETot = sqrt(GaussIntegralSE1^2 + GaussIntegralSE2^2);
    GaussSE1 = FitDeltas(1:6,:);
    GaussSE2 = FitDeltas(7:12,:);
    offsetSE = FitDeltas(end,:);
    centroid_mean = nanmean(centroidPositionVector);
    centroid_se = nanstd(centroidPositionVector);
    % combine for simplicity
    GaussIntegralVector = [GaussIntegral1 GaussIntegral2 GaussIntegralTot];
    GaussIntegralSEVector = [GaussIntegralSE1 GaussIntegralSE2 GaussIntegralSETot];
end