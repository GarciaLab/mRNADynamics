function  fitInfo = fit3DGaussian2spotV1(snip3D,PixelSize,zStep,varargin)
      
    % INPUT ARGUMENTS:
    % snip3D: 3D array containing spot to fit. Should contain only one spot
    % PixelSize: size of pixels (in nm)
    % zStep: pixels dimesnions along z direction (in nm)
    
    % RETURNS:
    %        fitInfo: data structure containing key fit results
    %
    % Parameter identity is as follows: 
    %        (1) amplitude of gaussian (spot 1)    
    %        (2-3) xy,and z sigma values (both spots) 
    %        (4-6) y,x,and z center positions (spot 1)         
    %        (7) amplitude of gaussian (spot 2)
    %        (8-10) y,x,and z center positions (spot 2)        
    %        (11-14) inferred background gradient            
    
    %% %%%%%% initialize inference params and define bounds %%%%%%%%%%%%%%%
    fitInfo = struct;
    fitInfo.fixedSpotDims = ~isempty(varargin);
    if ~isempty(varargin)
        spotDims = varargin{1};
    end
        
    % initial ballbark estimate for spot size
    % NL: these numbers reflect averages from a sample population of spots
    % In future, could be worthwhile performing some kind of hierarchical
    % fit to estimate population-wide PSFs
    sigmaXY_guess = 200/PixelSize;
    sigmaXY_guess_se = 40/PixelSize;
    sigmaZ_guess = 450/zStep;
    sigmaZ_guess_se = 175/zStep;
    
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
        mean(snip3D(:)/2),0,0,0]; % background fluorescence 
    
    % initialize upper and lower parameter bounds
    fitInfo.upperBoundVector = [...
        10*max(snip3D(:)), ... % Spot 1 amplitude
        sigmaXY_guess+sigmaXY_guess_se, sigmaZ_guess+sigmaZ_guess_se,... % Spot dimensions
        yDim+.5, xDim+.5, zDim+.5,... % Spot 1 position
        10*max(snip3D(:)),... % Spot 2 amplitude
        yDim+.5, xDim+.5, zDim+.5,... % Spot 2 position
        max(snip3D(:)), mean(snip3D(:)/2)/yDim, mean(snip3D(:)/2)/xDim, mean(snip3D(:)/2)/zDim]; % background fluorescence 
      
    fitInfo.lowerBoundVector = [...
        0, ... % Spot 1 amplitude
        sigmaXY_guess-sigmaXY_guess_se, sigmaZ_guess-sigmaZ_guess_se,... % Spot dimensions
        .5, .5, .5,... % Spot 1 position
        0,... % Spot 2 amplitude
        .5, .5, .5,... % Spot 2 position
        0, -mean(snip3D(:)/2)/yDim, -mean(snip3D(:)/2)/xDim, -mean(snip3D(:)/2)/zDim]; % background fluorescence    
    
       
    % define objective function
    fitInfo.dimensionVector = [yDim, xDim, zDim];    
    % define grid ref arrays
    [mesh_y, mesh_x, mesh_z] = meshgrid(1:yDim, 1:xDim, 1:zDim); 
    
    % define helper function to generate background fluo snip
    makeOffsetSnip = @(params) params(11) + params(12)*mesh_y + params(13)*mesh_x + params(14)*mesh_z;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% perform fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fitInfo.spot1ParamIndices = 1:6;
    fitInfo.spot2ParamIndices = [7 2:3 8:10];
    if ~fitInfo.fixedSpotDims 
        include_vec = true(1,14);        
        
        double3DObjective = @(params) simulate3DGaussSymmetric(fitInfo.dimensionVector, params(fitInfo.spot1ParamIndices))...
                                    + simulate3DGaussSymmetric(fitInfo.dimensionVector, params(fitInfo.spot2ParamIndices)) ...
                                    + makeOffsetSnip(params) - double(snip3D);           
    else
        include_vec = ~ismember(1:14,2:3);        
        
        double3DObjective = @(params) simulate3DGaussSymmetric(fitInfo.dimensionVector, [params(1) spotDims params(2:4)])...
                                    + simulate3DGaussSymmetric(fitInfo.dimensionVector, [params(5) spotDims params(6:8)]) ...
                                    + makeOffsetSnip([0 0 params]) - double(snip3D); 
    end
    
    
    % attempt to fit
    options.Display = 'off';
    [GaussFit, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(double3DObjective,...
        fitInfo.initial_parameters(include_vec),fitInfo.lowerBoundVector(include_vec),fitInfo.upperBoundVector(include_vec),options);
    
    % store parameters
    if fitInfo.fixedSpotDims         
        GaussFit = [GaussFit(1) spotDims GaussFit(2:end)];
    end
    fitInfo.RawFitParams = GaussFit;    
    fitInfo.Gauss1Params = GaussFit(fitInfo.spot1ParamIndices);
    fitInfo.Gauss2Params = GaussFit(fitInfo.spot2ParamIndices);    
  
    %% %%%%%%%%%%%%%%%%%%%%%% estimate uncertainty %%%%%%%%%%%%%%%%%%%%%%%%
    
    % estimate error in integral calculations numeriucally...this is faster
    % than doing it symbolically in matlab
    FitCI = nlparci(GaussFit(include_vec), residual, 'jacobian', jacobian);      
    FitDeltas = diff(FitCI') / 2 / 1.96;
    if fitInfo.fixedSpotDims       
        FitDeltas = [FitDeltas(1) 0 0 FitDeltas(2:end)];
    end
    fitInfo.RawFitSE = FitDeltas;
    gaussIntegral = @(params) params(1).*(2*pi).^1.5 .* params(2).^2.*params(3);  
    
    gaussIntegralError = @(params,paramsSE) (2*pi).^1.5 * sqrt(paramsSE(1)^2*(params(2)^2*params(3))^2 + ...
                            paramsSE(2)^2*(2*params(1)*params(2)*params(3))^2 + paramsSE(3)^2*(params(1)*params(2)^2)^2);                    
    
    % extract values to report
    fitInfo.GaussIntegral1 = gaussIntegral(GaussFit(fitInfo.spot1ParamIndices));
    fitInfo.GaussIntegralSE1 = gaussIntegralError(GaussFit(fitInfo.spot1ParamIndices),FitDeltas(fitInfo.spot1ParamIndices));
    fitInfo.GaussIntegral2 = gaussIntegral(GaussFit(fitInfo.spot2ParamIndices));
    fitInfo.GaussIntegralSE2 = gaussIntegralError(GaussFit(fitInfo.spot2ParamIndices),FitDeltas(fitInfo.spot2ParamIndices));
    fitInfo.GaussIntegralTot = fitInfo.GaussIntegral1 + fitInfo.GaussIntegral2;
    fitInfo.GaussIntegralSETot = sqrt(fitInfo.GaussIntegralSE1^2 + fitInfo.GaussIntegralSE2^2); % NL: it's likely that these two quantities are correlated, so addinin in quadrature is dubious
  
    fitInfo.Spot1Pos = GaussFit(4:6);
    fitInfo.Spot1PosSE = FitDeltas(4:6);
    fitInfo.Spot2Pos = GaussFit(8:10);
    fitInfo.Spot2PosSE = FitDeltas(8:10);
    fitInfo.SpotCentroid = (fitInfo.Spot1Pos*fitInfo.GaussIntegral1 + fitInfo.Spot2Pos*fitInfo.GaussIntegral2)/fitInfo.GaussIntegralTot;
    
    % estimate error in centroid position (also invokes dubious
    % independence assumption)
    int_vec = [fitInfo.GaussIntegral1 fitInfo.GaussIntegral2];
    se_vec = [fitInfo.GaussIntegralSE1 fitInfo.GaussIntegralSE2];
    fitInfo.SpotCentroidSE = sqrt(se_vec(1)^2.*(int_vec(2).*(fitInfo.Spot1Pos-fitInfo.Spot2Pos)./sum(int_vec).^2).^2 + ...
                                  fitInfo.Spot1PosSE.^2 * (int_vec(1)/sum(int_vec)).^2 + ...           
                                  se_vec(2)^2.*(int_vec(1).*(-fitInfo.Spot1Pos+fitInfo.Spot2Pos)./sum(int_vec).^2).^2 + ...
                                  fitInfo.Spot1PosSE.^2 * (int_vec(2)/sum(int_vec)).^2 ...
                                  );
                                
    fitInfo.offset = GaussFit(11)+GaussFit(12)*fitInfo.SpotCentroid(1)+GaussFit(13)*fitInfo.SpotCentroid(2)+GaussFit(14)*fitInfo.SpotCentroid(3);
    % calculate raw integral
    intMask = simulate3DGaussSymmetric(fitInfo.dimensionVector, [1 fitInfo.sigmaXY_int fitInfo.sigmaZ_int fitInfo.SpotCentroid]);    
    intMask = intMask >= exp(-9/2); % mask out everything > 3 sigma  
    offsetSnipFit = makeOffsetSnip(GaussFit);
%     offsetSnipVar = FitDeltas(11)^2 + FitDeltas(12)^2*mesh_y.^2 + FitDeltas(13)^2*mesh_x.^2 + FitDeltas(14)^2*mesh_z.^2;
    fitInfo.GaussIntegralRaw = (sum(intMask(:).*double(snip3D(:))) - sum(offsetSnipFit(:).*intMask(:)));
%     fitInfo.GaussIntegralRawSE = sqrt(sum(offsetSnipVar(:).*intMask(:)));
end