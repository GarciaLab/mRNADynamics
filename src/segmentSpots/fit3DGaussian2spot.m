function [GaussParams1, GaussParams2, offset, GaussIntegralVector, centroid_mean, GaussSE1, GaussSE2, offsetSE,...
    GaussIntegralSEVector, centroid_se]...
        = fit3DGaussian2spot(snip3D,PixelSize,zStep,varargin)
      
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
    %        (2-3) xy,and z sigma values (both spots) 
    %        (4-6) y,x,and z center positions (spot 1)         
    %        (7) amplitude of gaussian (spot 2)
    %        (8-10) y,x,and z center positions (spot 2)        
    %        (11) inferred background offset            
    
    %% %%%%%% initialize inference params and define bounds %%%%%%%%%%%%%%%
    fitInfo = struct;
    
    % initial ballbark estimate for spot size
    sigmaXY_guess = 200/PixelSize;
    sigmaZ_guess = 300/zStep;
    
    % define initial parameters
    xDim = size(snip3D,1);
    yDim = size(snip3D,2);
    zDim = size(snip3D,3);
    
    % initialize parameters
    fitInfo.initial_parameters =[...
        max(snip3D(:)), ... % Spot 1 amplitude
        sigmaXY_guess, sigmaZ_guess,... % Spot dimensions
        yDim/2-sigmaXY_guess,xDim/2-sigmaXY_guess, zDim/2-sigmaZ_guess,... % Spot 1 position
        max(snip3D(:)),... % Spot 2 amplitude
        yDim/2+sigmaXY_guess,xDim/2+sigmaXY_guess, zDim/2+sigmaZ_guess,... % Spot 2 position
        mean(snip3D(:)/2)]; % background fluorescence 
    
    % initialize upper and lower parameter bounds
    fitInfo.upperBoundVector = [...
        10*max(snip3D(:)), ... % Spot 1 amplitude
        5*sigmaXY_guess, 5*sigmaZ_guess,... % Spot dimensions
        yDim, xDim, zDim,... % Spot 1 position
        10*max(snip3D(:)),... % Spot 2 amplitude
        yDim, xDim, zDim,... % Spot 2 position
        max(snip3D(:))]; % background fluorescence 
      
    fitInfo.lowerBoundVector = [...
        0, ... % Spot 1 amplitude
        .25*sigmaXY_guess, .25*sigmaZ_guess,... % Spot dimensions
        1, 1, 1,... % Spot 1 position
        0,... % Spot 2 amplitude
        1, 1, 1,... % Spot 2 position
        0]; % background fluorescence    
    
       
    % define objective function
    fitInfo.dimensionVector = [yDim, xDim, zDim];    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% perform fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fitInfo.spot1ParamIndices = 1:6;
    fitInfo.spot2ParamIndices = [7 2:3 8:10];
    double3DObjective = @(params) simulate3DGaussSymmetric(fitInfo.dimensionVector, params(fitInfo.spot1ParamIndices))...
        + simulate3DGaussSymmetric(fitInfo.dimensionVector, params(fitInfo.spot2ParamIndices)) ...
        + params(end) - double(snip3D);           
    
    % attempt to fit
    options.Display = 'off';
    [GaussFit, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(double3DObjective,...
        fitInfo.initial_parameters,fitInfo.lowerBoundVector,fitInfo.upperBoundVector,options);
      
    % store parameters
    fitInfo.Gauss1Params = GaussFit(fitInfo.spot1ParamIndices);
    fitInfo.Gauss2Params = GaussFit(fitInfo.spot2ParamIndices);
    fitInfo.offset = GaussFit(end);
    
    %% %%%%%%%%%%%%%%%%%%%%%% estimate uncertainty %%%%%%%%%%%%%%%%%%%%%%%%
    
    % estimate error in integral calculations numeriucally
    FitCI = nlparci(GaussFit, residual, 'jacobian', jacobian);
    FitCI(FitCI<0) = 0; % Meh?
    FitDeltas = diff(FitCI')' / 2 / 1.96;
    paramValues = normrnd(0,1,size(FitDeltas,1),100).*FitDeltas + GaussFit';
    paramValues(paramValues<0) = realmin;
    
    % de
    covarianceMatrix = @(params) [params(2),    0,    0;
                                  0,    params(2),    0;
                                  0,    0,     params(3)];    
                                
    gaussIntegral = @(params) params(1)*(2*pi)^1.5 * sqrt(det(covarianceMatrix(params)));    
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