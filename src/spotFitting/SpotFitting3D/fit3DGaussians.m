function  fitInfo = fit3DGaussians(snip3D, PixelSize, zStep, nSpots, givenParams)

% DESCRIPTION
%
%
% INPUT ARGUMENTS
% snip3D: 3D array containing spot to fit. Should contain only one spot
% PixelSize: pixel size in later (xy) dimension (in nm)
% zStep: pixels size in axial (z) direction (in nm)
% nSpots: number of gaussians to fit. Should be 2 for MS2 spots (to
%         account for sister chromatids) and 1 for transcription factor
%         clusters
% givenParams:
%
% OPTIONS
% N/A
%
% OUTPUT
% fitInfo: Data structure containing key fit results.
%          Parameter identity is as follows: 
%               (1) amplitude of gaussian (spot 1)    
%               (2-3) xy, and z sigma values (both spots) 
%               (4-6) y,x,and z center positions (spot 1)         
%               (7) amplitude of gaussian (spot 2)
%               (8-10) y,x,and z center positions (spot 2)        
%               (11-17) inferred background gradient
%
% Author (contact): Nicholas Lammers (nlammers@berkeley.edu)
% Created: 2019/2020ish
% Last Updated:

%% %%%%%% initialize inference params and define bounds %%%%%%%%%%%%%%%
fitInfo = struct;
fitInfo.twoSpotFlag = (nSpots==2);    
fitInfo.fitFlag = isempty(givenParams);

% NL: these numbers reflect averages from a sample population of spots
% In future, could be worthwhile performing some kind of hierarchical
% fit to estimate population-wide PSFs
sigmaXY_guess = 200/PixelSize;
sigmaXY_guess_se = 100/PixelSize;
sigmaZ_guess = 450/zStep;
sigmaZ_guess_se = 225/zStep;

% set size of neighborhood to integrate for raw integral
fitInfo.sigmaXY_int = 230 / PixelSize;
fitInfo.sigmaZ_int = 620 / zStep;

% define initial parameters
xDim = size(snip3D,1);
yDim = size(snip3D,2);
zDim = size(snip3D,3);      

% initialize upper and lower parameter bounds
xyPosBound = 0.5;
zPosBound = 0.5;
bkgdFluorBound = mean(snip3D(:)/2); % MT: why are are defining it this way?

maxPixIntensity = max(snip3D(:));
spotAmpUpperBound = 4*maxPixIntensity;

fitInfo.upperBoundVector = [...
    spotAmpUpperBound, ... % (1) Spot 1 amplitude
    sigmaXY_guess+sigmaXY_guess_se, sigmaZ_guess+sigmaZ_guess_se,... % (2-3) Spot dimensions
    yDim+xyPosBound, xDim+xyPosBound, zDim+zPosBound,... % (4-6) Spot 1 position
    spotAmpUpperBound,... % (7) Spot 2 amplitude
    yDim+xyPosBound, xDim+xyPosBound, zDim+zPosBound,... % (8-10) Spot 2 position
    maxPixIntensity, ...
    bkgdFluorBound/yDim, bkgdFluorBound/xDim, bkgdFluorBound/zDim,...        
    bkgdFluorBound/yDim, bkgdFluorBound/xDim, bkgdFluorBound/zDim]; % (11-17) background fluorescence 

minPixIntensity = 0; 
spotAmpLowerBound = bkgdFluorBound;

fitInfo.lowerBoundVector = [...
    spotAmpLowerBound, ... % (1) Spot 1 amplitude
    sigmaXY_guess-sigmaXY_guess_se, sigmaZ_guess-sigmaZ_guess_se,... % (2-3) Spot dimensions
    xyPosBound, xyPosBound, zPosBound,... % (4-6) Spot 1 position
    spotAmpLowerBound,... % (7) Spot 2 amplitude
    xyPosBound, xyPosBound, zPosBound,... % (8-10) Spot 2 position
    minPixIntensity,...
    -bkgdFluorBound/yDim, -bkgdFluorBound/xDim, -bkgdFluorBound/zDim,...
    -bkgdFluorBound/yDim, -bkgdFluorBound/xDim, -bkgdFluorBound/zDim]; % (11-17) background fluorescence    


% Initialize Parameters
% If we expect only 1 spot, initialize spot1 position in the center of the
% snip, but if we expect 2 spots, initialize position offset from center
fitInfo.initial_parameters =[...
    maxPixIntensity, ... % (1) Spot 1 amplitude
    sigmaXY_guess, sigmaZ_guess,... % (2-3) Spot dimensions
    (yDim/2 - sigmaXY_guess*fitInfo.twoSpotFlag),...
    (xDim/2 - sigmaXY_guess*fitInfo.twoSpotFlag),...
    (zDim/2 - sigmaZ_guess*fitInfo.twoSpotFlag),...% (4-6) Spot 1 position
    maxPixIntensity,... % (7) Spot 2 amplitude
    yDim/2+sigmaXY_guess, ...
    xDim/2+sigmaXY_guess, ...
    zDim/2+sigmaZ_guess,... % (8-10) Spot 2 position
    bkgdFluorBound,...
    0,0,0,...
    0,0,0]; % (11-17) background fluorescence 


% define objective function
fitInfo.dimensionVector = [yDim, xDim, zDim];  

% define grid ref arrays
[mesh_y, mesh_x, mesh_z] = meshgrid(1:yDim, 1:xDim, 1:zDim); 

% define helper function to generate background fluo snip
if fitInfo.twoSpotFlag
    makeOffsetSnip = @(params) params(11) + ...
                               params(12)*mesh_y + params(13)*mesh_x + params(14)*mesh_z + ...
                               params(15)*mesh_y.^2 + params(16)*mesh_x.^2 + params(17)*mesh_z.^2;

    makeOffsetSnipSE = @(params)  sqrt(params(11)^2 );% + ... NL: gradient parameters are generally small and somehwat volatile, so they lead to inflated error estimates                                                                            
%                                            params(12)^2*mesh_y.^2 + params(13)^2*mesh_x.^2 + params(14)^2*mesh_z.^2 + ...
%                                            params(15)^2*mesh_y.^4 + params(16)^2*mesh_x.^4 + params(17)^2*mesh_z.^4);
    include_vec = true(1,17); 
else
    makeOffsetSnip = @(params) params(7) + ...
                               params(8)*mesh_y + params(9)*mesh_x + params(10)*mesh_z + ...
                               params(11)*mesh_y.^2 + params(12)*mesh_x.^2 + params(13)*mesh_z.^2;

    makeOffsetSnipSE = @(params)  params(7)^2 ;
    
    % exclude spot2 parameters if we're only fitting one spot
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
    fitInfo.GaussIntegralSETot = sqrt(fitInfo.GaussIntegralSE1^2 + fitInfo.GaussIntegralSE2^2); % NL: it's likely that these two quantities are correlated, so addinin in quadrature is dubious
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
offsetSnipSE = makeOffsetSnipSE(FitDeltas);
fitInfo.GaussIntegralRaw = (sum(intMask(:).*double(snip3D(:)) - offsetSnipFit(:).*intMask(:)));
fitInfo.GaussIntegralRawSE = sum(offsetSnipSE*intMask(:));
if fitInfo.twoSpotFlag
    fitInfo.offset = GaussFit(11) + ...
                     GaussFit(12:14)*fitInfo.SpotCentroid' + ...
                     GaussFit(15:17)*((fitInfo.SpotCentroid').^2);

else
    fitInfo.offset = GaussFit(7) + ...
                     GaussFit(8:10)*fitInfo.SpotCentroid' + ...
                     GaussFit(11:13)*((fitInfo.SpotCentroid').^2);
end