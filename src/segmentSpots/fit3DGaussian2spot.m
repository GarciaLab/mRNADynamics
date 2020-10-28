function [GaussParams1, GaussParams2, offset, GaussIntegralVector, spotCentroid, GaussSE1, GaussSE2, offsetSE,...
    GaussIntegralSEVector, spotCentroidSE]...
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
    sigmaXY_guess = 100/PixelSize;
    sigmaZ_guess = 300/zStep;
        
    % set size of nerighborhood to integrate for raw integral
    fitInfo.sigmaXY_int = 230 / PixelSize;
    fitInfo.sigmaZ_int = 620 / zStep;
    
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
    
    % estimate error in integral calculations numeriucally...this is faster
    % than doing it symbolically in matlab
    FitCI = nlparci(GaussFit, residual, 'jacobian', jacobian);    
    FitDeltas = diff(FitCI')' / 2 / 1.96;
    nSamples = 100;
    paramValues = normrnd(0,1,size(FitDeltas,1),nSamples).*FitDeltas + GaussFit';
    paramValues(paramValues<0) = realmin;
    
    % define helper functions for integral calculation
    covarianceMatrix = @(params) [params(2)^2,    0,    0;
                                  0,      params(2)^2,    0;
                                  0,    0,        params(3)^2];    
                                
    gaussIntegral = @(params) params(1)*(2*pi)^1.5 * sqrt(det(covarianceMatrix(params)));    
    gaussIntegralSpot1 = NaN(1,nSamples);
    gaussIntegralSpot2 = NaN(1,nSamples);    
    gaussIntegralTotal = NaN(1,nSamples);    
    positionVectorSpot1 = NaN(nSamples,3); 
    positionVectorSpot2 = NaN(nSamples,3); 
    centroidPositionVector = NaN(nSamples,3); 
    
    for i = 1:size(paramValues,2)
        gaussIntegralSpot1(i) = gaussIntegral(paramValues(fitInfo.spot1ParamIndices,i));
        gaussIntegralSpot2(i) = gaussIntegral(paramValues(fitInfo.spot2ParamIndices,i));  
        
        gaussIntegralTotal(i) = gaussIntegralSpot1(i) + gaussIntegralSpot2(i);
        
        positionVectorSpot1(i,:) = paramValues(fitInfo.spot1ParamIndices(4:6),i);
        positionVectorSpot2(i,:) = paramValues(fitInfo.spot2ParamIndices(4:6),i);
        
        centroidPositionVector(i,:) = (gaussIntegralSpot1(i)*positionVectorSpot1(i,:) + ...
                                       gaussIntegralSpot2(i)*positionVectorSpot2(i,:))/...
                                       gaussIntegralTotal(i);
    end
    
    % extract values to report
    fitInfo.GaussIntegral1 = mean(gaussIntegralSpot1);
    fitInfo.GaussIntegralSE1 = std(gaussIntegralSpot1);
    fitInfo.GaussIntegral2 = mean(gaussIntegralSpot2);
    fitInfo.GaussIntegralSE2 = std(gaussIntegralSpot2);
    fitInfo.GaussIntegralTot = mean(gaussIntegralTotal);
    fitInfo.GaussIntegralSETot = std(gaussIntegralTotal);

    fitInfo.Spot1Pos = mean(positionVectorSpot1);
    fitInfo.Spot1PosSE = std(positionVectorSpot1);
    fitInfo.Spot2Pos = mean(positionVectorSpot2);
    fitInfo.Spot2PosSE = std(positionVectorSpot2);
    fitInfo.spotCentroid = mean(centroidPositionVector);
    fitInfo.spotCentroidSE = std(centroidPositionVector);
    
    % combine for simplicity
    GaussIntegralVector = [GaussIntegral1 GaussIntegral2 GaussIntegralTot];
    GaussIntegralSEVector = [GaussIntegralSE1 GaussIntegralSE2 GaussIntegralSETot];
    
    % calculate raw integral
    intMask = simulate3DGaussSymmetric(fitInfo.dimensionVector, [1 fitInfo.sigmaXY_int fitInfo.sigmaZ_int fitInfo.spotCentroid]);
    intMask = intMask / sum(intMask(:));
    [intSorted,~] = sort(intMask(:),'descend');
    cutoffVal = intSorted(find(cumsum(intSorted)>0.997,1));
    intMask = intMask >= cutoffVal;    
    fitInfo.GaussIntegralRaw = (sum(intMask(:).*double(snip3D(:))) - sum(GaussFit(end).*intMask(:)));
end