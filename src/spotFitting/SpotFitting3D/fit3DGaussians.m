function  fitInfo = fit3DGaussians(snip3D, PixelSize, zStep, spotDims, nSpots, givenParams)
      
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
    %        (11-17) inferred background gradient            
    
    %% %%%%%% initialize inference params and define bounds %%%%%%%%%%%%%%%
    fitInfo = struct;    
    fitInfo.spotInfoFlag = ~isempty(spotDims);
    fitInfo.spotDims = spotDims;
    fitInfo.twoSpotFlag = nSpots==2;    
    fitInfo.fitFlag = isempty(givenParams);
    
    % initial ballbark estimate for spot size
    if fitInfo.spotInfoFlag
        sigmaXY_guess = spotDims.sigmaXY;
        sigmaXY_guess_se = spotDims.sigmaXYSE;
        sigmaZ_guess = spotDims.sigmaZ;
        sigmaZ_guess_se = spotDims.sigmaZSE;
    else  
        % NL: these numbers reflect averages from a sample population of spots
        % In future, could be worthwhile performing some kind of hierarchical
        % fit to estimate population-wide PSFs
        sigmaXY_guess = 200/PixelSize;
        sigmaXY_guess_se = 100/PixelSize;
        sigmaZ_guess = 450/zStep;
        sigmaZ_guess_se = 225/zStep;
    end
    % set size of nerighborhood to integrate for raw integral
    fitInfo.sigmaXY_int = 230 / PixelSize;
    fitInfo.sigmaZ_int = 620 / zStep;
    
    % define initial parameters
    xDim = size(snip3D,1);
    yDim = size(snip3D,2);
    zDim = size(snip3D,3);      
    
    % initialize upper and lower parameter bounds
    fitInfo.upperBoundVector = [...
        4*max(snip3D(:)), ... % Spot 1 amplitude
        sigmaXY_guess+sigmaXY_guess_se, sigmaZ_guess+sigmaZ_guess_se,... % Spot dimensions
        yDim+.5, xDim+.5, zDim+.5,... % Spot 1 position
        4*max(snip3D(:)),... % Spot 2 amplitude
        yDim+.5, xDim+.5, zDim+.5,... % Spot 2 position
        max(snip3D(:)), ...
        mean(snip3D(:)/2)/yDim, mean(snip3D(:)/2)/xDim, mean(snip3D(:)/2)/zDim,...        
        mean(snip3D(:)/2)/yDim, mean(snip3D(:)/2)/xDim, mean(snip3D(:)/2)/zDim]; % background fluorescence 
      
    fitInfo.lowerBoundVector = [...
        mean(snip3D(:)/2), ... % Spot 1 amplitude
        sigmaXY_guess-sigmaXY_guess_se, sigmaZ_guess-sigmaZ_guess_se,... % Spot dimensions
        .5, .5, .5,... % Spot 1 position
        mean(snip3D(:)/2),... % Spot 2 amplitude
        .5, .5, .5,... % Spot 2 position
        0,...
        -mean(snip3D(:)/2)/yDim, -mean(snip3D(:)/2)/xDim, -mean(snip3D(:)/2)/zDim,...
        -mean(snip3D(:)/2)/yDim, -mean(snip3D(:)/2)/xDim, -mean(snip3D(:)/2)/zDim]; % background fluorescence    
    

    % initialize parameters
    fitInfo.initial_parameters =[...
        max(snip3D(:)), ... % Spot 1 amplitude
        sigmaXY_guess, sigmaZ_guess,... % Spot dimensions
        yDim/2-sigmaXY_guess*fitInfo.twoSpotFlag,xDim/2-sigmaXY_guess*fitInfo.twoSpotFlag, zDim/2-sigmaZ_guess*fitInfo.twoSpotFlag,... % Spot 1 position
        max(snip3D(:)),... % Spot 2 amplitude
        yDim/2+sigmaXY_guess,xDim/2+sigmaXY_guess, zDim/2+sigmaZ_guess,... % Spot 2 position
        mean(snip3D(:)/2),...
        0,0,0,...
        0,0,0]; % background fluorescence 
      
    
    % define objective function
    fitInfo.dimensionVector = [yDim, xDim, zDim];  
    
    % define grid ref arrays
    [mesh_y, mesh_x, mesh_z] = meshgrid(1:yDim, 1:xDim, 1:zDim); 
    
    % define helper function to generate background fluo snip
    if fitInfo.twoSpotFlag
        makeOffsetSnip = @(params) params(11) + ...
                                   params(12)*mesh_y + params(13)*mesh_x + params(14)*mesh_z + ...
                                   params(15)*mesh_y.^2 + params(16)*mesh_x.^2 + params(17)*mesh_z.^2;
        include_vec = true(1,17); 
    else
        makeOffsetSnip = @(params) params(7) + ...
                                   params(8)*mesh_y + params(9)*mesh_x + params(10)*mesh_z + ...
                                   params(11)*mesh_y.^2 + params(12)*mesh_x.^2 + params(13)*mesh_z.^2;
        include_vec = ~ismember(1:17,7:10); 
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% perform fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fitInfo.spot1ParamIndices = 1:6;
    fitInfo.spot2ParamIndices = [7 2:3 8:10];

           
    if fitInfo.twoSpotFlag        
        spot3DObjective = @(params) simulate3DGaussSymmetric(fitInfo.dimensionVector, params(fitInfo.spot1ParamIndices))...
                                    + simulate3DGaussSymmetric(fitInfo.dimensionVector, params(fitInfo.spot2ParamIndices)) ...
                                    + makeOffsetSnip(params) - double(snip3D);               
    else 
        spot3DObjective = @(params) simulate3DGaussSymmetric(fitInfo.dimensionVector, params(fitInfo.spot1ParamIndices))...
                                    + makeOffsetSnip(params) - double(snip3D);     
    end 
    
    options.Display = 'off';   
    
    % update initialization and bound fields
    fitInfo.initial_parameters = fitInfo.initial_parameters(include_vec);
    fitInfo.upperBoundVector = fitInfo.upperBoundVector(include_vec);
    fitInfo.lowerBoundVector = fitInfo.lowerBoundVector(include_vec);
    
    % attempt to fit
    if fitInfo.fitFlag
        [GaussFit, fitInfo.resnorm, residual, fitInfo.exitflag, ~,~,jacobian] = lsqnonlin(spot3DObjective,...
            fitInfo.initial_parameters,fitInfo.lowerBoundVector,fitInfo.upperBoundVector,options);
    else
        options.MaxIterations = 1;
        [~, fitInfo.resnorm, residual, fitInfo.exitflag, ~,~,jacobian] = lsqnonlin(spot3DObjective,...
            givenParams(include_vec),fitInfo.lowerBoundVector,fitInfo.upperBoundVector,options);
        GaussFit = givenParams(include_vec);
    end
        
    % require that spot1 is always closer to the origin
    orderIndices = 1:length(GaussFit);
    if fitInfo.twoSpotFlag  
        spot1Score = sum(GaussFit(4:6).^2);
        spot2Score = sum(GaussFit(8:10).^2);        
        if spot1Score > spot2Score        
            orderIndices(fitInfo.spot1ParamIndices) = fitInfo.spot2ParamIndices;
            orderIndices(fitInfo.spot2ParamIndices) = fitInfo.spot1ParamIndices;
        end           
    end
    %% %%%%%%%%%%%%%%%%%%%%%% estimate uncertainty %%%%%%%%%%%%%%%%%%%%%%%%
    
    % estimate error in integral calculations numeriucally...this is faster
    % than doing it symbolically in matlab    
    FitCI = nlparci(GaussFit, residual, 'jacobian', jacobian);      
    FitDeltas = diff(FitCI') / 2 / 1.96;
    
    % Reorder if necessary
    GaussFit = GaussFit(orderIndices);
    FitDeltas = FitDeltas(orderIndices);
    
    % store parameters
    fitInfo.RawFitParams = GaussFit;    
    fitInfo.Gauss1Params = GaussFit(fitInfo.spot1ParamIndices);
    if fitInfo.twoSpotFlag
        fitInfo.Gauss2Params = GaussFit(fitInfo.spot2ParamIndices);    
    else
        fitInfo.Gauss2Params = NaN(size(fitInfo.spot2ParamIndices));
    end
    
    fitInfo.RawFitSE = FitDeltas;
    gaussIntegral = @(params) params(1).*(2*pi).^1.5 .* params(2).^2.*params(3);  
    
    gaussIntegralError = @(params,paramsSE) (2*pi).^1.5 * sqrt(paramsSE(1)^2*(params(2)^2*params(3))^2 + ...
                            paramsSE(2)^2*(2*params(1)*params(2)*params(3))^2 + paramsSE(3)^2*(params(1)*params(2)^2)^2);                    
    
    % extract values to report
    fitInfo.GaussIntegral1 = gaussIntegral(GaussFit(fitInfo.spot1ParamIndices));
    fitInfo.GaussIntegralSE1 = gaussIntegralError(GaussFit(fitInfo.spot1ParamIndices),FitDeltas(fitInfo.spot1ParamIndices));
    fitInfo.Spot1Pos = GaussFit(4:6);
    fitInfo.Spot1PosSE = FitDeltas(4:6);
    
    if fitInfo.twoSpotFlag
        fitInfo.GaussIntegral2 = gaussIntegral(GaussFit(fitInfo.spot2ParamIndices));
        fitInfo.GaussIntegralSE2 = gaussIntegralError(GaussFit(fitInfo.spot2ParamIndices),FitDeltas(fitInfo.spot2ParamIndices));
        fitInfo.GaussIntegralTot = fitInfo.GaussIntegral1 + fitInfo.GaussIntegral2;
        fitInfo.GaussIntegralSETot = NaN;%sqrt(fitInfo.GaussIntegralSE1^2 + fitInfo.GaussIntegralSE2^2); % NL: it's likely that these two quantities are correlated, so addinin in quadrature is dubious
        fitInfo.Spot2Pos = GaussFit(8:10);
        fitInfo.Spot2PosSE = FitDeltas(8:10);
        fitInfo.SpotCentroid = (fitInfo.Spot1Pos*fitInfo.GaussIntegral1^2 + fitInfo.Spot2Pos*fitInfo.GaussIntegral2^2)/...
                               (fitInfo.GaussIntegral1^2 + fitInfo.GaussIntegral2^2);%fitInfo.GaussIntegralTot;
        
        % estimate error in centroid position (also invokes dubious
        % independence assumption)
        int_vec = [fitInfo.GaussIntegral1 fitInfo.GaussIntegral2];
        se_vec = [fitInfo.GaussIntegralSE1 fitInfo.GaussIntegralSE2];
        fitInfo.SpotCentroidSE = sqrt(se_vec(1)^2.*(int_vec(2).*(fitInfo.Spot1Pos-fitInfo.Spot2Pos)./sum(int_vec).^2).^2 + ...
                                      fitInfo.Spot1PosSE.^2 * (int_vec(1)/sum(int_vec)).^2 + ...           
                                      se_vec(2)^2.*(int_vec(1).*(-fitInfo.Spot1Pos+fitInfo.Spot2Pos)./sum(int_vec).^2).^2 + ...
                                      fitInfo.Spot1PosSE.^2 * (int_vec(2)/sum(int_vec)).^2 ...
                                      );
    else
        fitInfo.GaussIntegral2 = NaN;
        fitInfo.GaussIntegralSE2 = NaN;
        fitInfo.GaussIntegralTot = fitInfo.GaussIntegral1;
        fitInfo.GaussIntegralSETot = fitInfo.GaussIntegralSE1;
        fitInfo.Spot2Pos = NaN(1,3);
        fitInfo.Spot2PosSE = NaN(1,3);
        fitInfo.SpotCentroid = fitInfo.Spot1Pos;
        fitInfo.SpotCentroidSE = fitInfo.Spot1PosSE;
    end              
        
    
    % calculate raw integral
    cutoffVal = exp(-6.25/2); %2.5 sigma
    if fitInfo.twoSpotFlag
        intMask1 = simulate3DGaussSymmetric(fitInfo.dimensionVector, [1 fitInfo.sigmaXY_int fitInfo.sigmaZ_int fitInfo.Spot1Pos]);    
        intMask2 = simulate3DGaussSymmetric(fitInfo.dimensionVector, [1 fitInfo.sigmaXY_int fitInfo.sigmaZ_int fitInfo.Spot2Pos]);    
        intMask = intMask1 >= cutoffVal| intMask2 >= cutoffVal;  
    else
        intMask = simulate3DGaussSymmetric(fitInfo.dimensionVector, [1 fitInfo.sigmaXY_int fitInfo.sigmaZ_int fitInfo.Spot1Pos]);    
        intMask = intMask >= cutoffVal; 
    end
    
    offsetSnipFit = makeOffsetSnip(GaussFit);
    fitInfo.GaussIntegralRaw = (sum(intMask(:).*double(snip3D(:)) - offsetSnipFit(:).*intMask(:)));

    if fitInfo.twoSpotFlag
        fitInfo.offset = GaussFit(11) + ...
                         GaussFit(12:14)*fitInfo.SpotCentroid' + ...
                         GaussFit(15:17)*((fitInfo.SpotCentroid').^2);

    else
        fitInfo.offset = GaussFit(7) + ...
                         GaussFit(8:10)*fitInfo.SpotCentroid' + ...
                         GaussFit(11:13)*((fitInfo.SpotCentroid').^2);
    end