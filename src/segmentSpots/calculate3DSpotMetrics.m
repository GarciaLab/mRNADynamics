function fitInfo = calculate3DSpotMetrics(fitInfo, snip3D)
        
    
    GaussFit = fitInfo.RawFitParams;
    
    gaussIntegral = @(params) params(1).*(2*pi).^1.5 .* params(2).^2.*params(3);  
       
    % extract values to report
    fitInfo.GaussIntegral1 = gaussIntegral(GaussFit(fitInfo.spot1ParamIndices));
   
    fitInfo.Spot1Pos = GaussFit(4:6);
   
    if fitInfo.twoSpotFlag
        fitInfo.GaussIntegral2 = gaussIntegral(GaussFit(fitInfo.spot2ParamIndices));
        
        fitInfo.GaussIntegralTot = fitInfo.GaussIntegral1 + fitInfo.GaussIntegral2;
        
        fitInfo.Spot2Pos = GaussFit(8:10);
        
        fitInfo.SpotCentroid = (fitInfo.Spot1Pos*fitInfo.GaussIntegral1 + fitInfo.Spot2Pos*fitInfo.GaussIntegral2)/fitInfo.GaussIntegralTot;
                
        
    else
        fitInfo.GaussIntegral2 = NaN;        
        fitInfo.GaussIntegralTot = fitInfo.GaussIntegral1;
        
        fitInfo.Spot2Pos = NaN(1,3);
        
        fitInfo.SpotCentroid = fitInfo.Spot1Pos;
        
    end              
        
    if fitInfo.twoSpotFlag
        fitInfo.offset = GaussFit(11)+GaussFit(12)*fitInfo.SpotCentroid(1)+GaussFit(13)*fitInfo.SpotCentroid(2)+GaussFit(14)*fitInfo.SpotCentroid(3);
    else
        fitInfo.offset = GaussFit(7)+GaussFit(8)*fitInfo.SpotCentroid(1)+GaussFit(9)*fitInfo.SpotCentroid(2)+GaussFit(10)*fitInfo.SpotCentroid(3);
    end
    
    % define reference arrays
    [mesh_y, mesh_x, mesh_z] = meshgrid(1:fitInfo.dimensionVector(1), 1:fitInfo.dimensionVector(2), 1:fitInfo.dimensionVector(3)); 
    
    % define helper function to generate background fluo snip
    if fitInfo.twoSpotFlag
        makeOffsetSnip = @(params) params(11) + params(12)*mesh_y + params(13)*mesh_x + params(14)*mesh_z;        
    else
        makeOffsetSnip = @(params) params(7) + params(8)*mesh_y + params(9)*mesh_x + params(10)*mesh_z;        
    end
    
    
    % calculate raw integral
    intMask = simulate3DGaussSymmetric(fitInfo.dimensionVector, [1 fitInfo.sigmaXY_int fitInfo.sigmaZ_int fitInfo.SpotCentroid]);    
    intMask = intMask >= exp(-9/2); % mask out everything > 3 sigma  
    offsetSnipFit = makeOffsetSnip(GaussFit);
%     offsetSnipVar = FitDeltas(11)^2 + FitDeltas(12)^2*mesh_y.^2 + FitDeltas(13)^2*mesh_x.^2 + FitDeltas(14)^2*mesh_z.^2;
    fitInfo.GaussIntegralRaw = (sum(intMask(:).*double(snip3D(:))) - sum(offsetSnipFit(:).*intMask(:)));