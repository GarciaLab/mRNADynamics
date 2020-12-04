clear
close all

% set some basic simulation parameters
snipSize = [29 29 12];
stackSize = [50 50 14];
nSpotsToSimulate = 100;

% load Spots data for reference set
Prefix = '2020-03-15-2xDl-Ven_hbP2P-mCh_Leica_Zoom2_7uW14uW_05';
liveExperiment = LiveExperiment(Prefix);
Spots = getSpots(liveExperiment);

%% store info from individual fits to 
paramArray = [];
iter = 1;
for s = 1:length(Spots)
    for f = 1:length(Spots(s).Fits)
        paramArray(iter,:) = Spots(s).Fits(f).SpotFitInfo3D.RawFitParams;
        iter = iter + 1;
    end
end

% generate ranges for the different parameter values
paramVecUb = prctile(paramArray,90);
paramVecLb = prctile(paramArray,10);
paramDelta = paramVecUb-paramVecLb;
% define helper functions
spot1ParamIndices = 1:6;
spot2ParamIndices = [7 2:3 8:10];

[mesh_y, mesh_x, mesh_z] = meshgrid(1:snipSize(1), 1:snipSize(2), 1:snipSize(3)); 

makeOffsetSnip = @(params) params(11) + params(12)*mesh_y + params(13)*mesh_x + params(14)*mesh_z;
        include_vec = true(1,14); 
    
spot3DFun = @(params) simulate3DGaussSymmetric(snipSize, params(spot1ParamIndices))...
                                    + simulate3DGaussSymmetric(snipSize, params(spot2ParamIndices)) ...
                                    + makeOffsetSnip(params);
  
                                  
% generate simulated stacks
simStackCell = cell(1,nSpotsToSimulate);
% store paramter values
trueParamArray = NaN(nSpotsToSimulate,length(paramDelta));
infParamArray = NaN(nSpotsToSimulate,length(paramDelta));

startParallelPool(8, 0, 1);
parfor n = 1:nSpotsToSimulate
    % draw parameter values
    trueParamArray(n,:) = rand(size(paramDelta)).*paramDelta + paramVecLb;
    % create simulated snip
    simStackCell{n} = spot3DFun(simParams);
    % infer parameters back
    fitInfo = fit3DGaussians(simStackCell{n},liveExperiment.pixelSize_nm,liveExperiment.zStep_um*1000,[],2);
    infParamArray(n,:) = fitInfo.RawFitParams;
end
    
  
        