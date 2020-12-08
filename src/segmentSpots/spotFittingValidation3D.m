clear
close all

% set some basic simulation parameters
snipSize = [29 29 12];
stackSize = [50 50 14];
nSpotsToSimulate = 5e2;
nBootstraps = 10;
% load Spots data for reference set
Prefix = '2020-03-15-2xDl-Ven_hbP2P-mCh_Leica_Zoom2_7uW14uW_05';
liveExperiment = LiveExperiment(Prefix);
FrameInfo = getFrameInfo(liveExperiment);
Spots = getSpots(liveExperiment);

% get pixel dimensions
PixelsSizeXY = FrameInfo(1).PixelSize;
PixelsSizeZ = FrameInfo(1).ZStep;

%% store info from individual fits to 
paramArray = [];
iter = 1;
for s = 1:length(Spots)
    for f = 1:length(Spots(s).Fits)
        paramArray(iter,:) = Spots(s).Fits(f).SpotFitInfo3D.RawFitParams;
        iter = iter + 1;
    end
end
% reparameterize spot 2 position wrpt spot 1
spotDistArray = paramArray(:,4:6)-paramArray(:,8:10);

% generate ranges for the different parameter values
paramVecUb = prctile(paramArray,90);
paramVecLb = prctile(paramArray,10);
paramDelta = paramVecUb-paramVecLb;

% define helper functions
spot1ParamIndices = 1:6;
spot2ParamIndices = [7 2:3 8:10];

[mesh_y, mesh_x, mesh_z] = meshgrid(1:snipSize(1), 1:snipSize(2), 1:snipSize(3)); 

% define function to make noisy snip background
makeOffsetSnip = @(params) normrnd(0,.5*params(11),size(mesh_z)) + params(11) + params(12)*mesh_y + params(13)*mesh_x + params(14)*mesh_z;        
    
spot3DFun = @(params) simulate3DGaussSymmetric(snipSize, params(spot1ParamIndices))...
                                    + simulate3DGaussSymmetric(snipSize, params(spot2ParamIndices)) ...
                                    + makeOffsetSnip(params);
  
%initialize structure to store validation results
validationInfo = struct;
validationInfo(nSpotsToSimulate) = struct;

% store paramter values

startParallelPool(24, 0, 1);

disp('Simulating spots and fitting...')
parfor n = 1:nSpotsToSimulate
    % draw parameter values
    trueParams = rand(size(paramDelta)).*paramDelta + paramVecLb;
    validationInfo(n).trueParams = trueParams;
    
    % create simulated snip
    simStack = spot3DFun(trueParams);
    validationInfo(n).simStack = simStack;
    
    % infer parameters back
    tempParams = NaN(nBootstraps,14);
    residVec = NaN(1,nBootstraps);
    for b = 1:nBootstraps
        fitInfoRaw = fit3DGaussians(validationInfo(n).simStack,liveExperiment.pixelSize_nm,liveExperiment.zStep_um*1000,[],2,b~=nBootstraps);
        validationInfo(n).fitInfo{b} = fitInfoRaw;
        tempParams(b,:) = fitInfoRaw.RawFitParams;
        residVec(b) = fitInfoRaw.resnorm;
    end        
    
    % botain true values for key spot parameters
%     trueInfo.RawFitParams = fitInfo.;
    trueInfo = struct;
    trueInfo.spot1ParamIndices = fitInfoRaw.spot1ParamIndices;
    trueInfo.spot2ParamIndices = fitInfoRaw.spot2ParamIndices;
    trueInfo.twoSpotFlag = fitInfoRaw.twoSpotFlag;
    trueInfo.dimensionVector = fitInfoRaw.dimensionVector;
    trueInfo.sigmaXY_int = fitInfoRaw.sigmaXY_int;
    trueInfo.sigmaZ_int = fitInfoRaw.sigmaZ_int;
    
    trueInfo.RawFitParams = trueParams;
    
    trueInfo = calculate3DSpotMetrics(trueInfo,simStack);

    % calculate average for results         
    fitInfoBoot = trueInfo;
    fitInfoBoot.RawFitParams = mean(tempParams);    
    % initialize structure for bootstrap mean
    fitInfoBoot = calculate3DSpotMetrics(fitInfoBoot,simStack);
    validationInfo(n).fitInfoBoot = fitInfoBoot;
    
    % calculate ml for results         
    fitInfoML = trueInfo;
    [~,mi] = min(residVec);
    fitInfoML.RawFitParams = tempParams(mi,:);    
    % initialize structure for bootstrap mean
    fitInfoML = calculate3DSpotMetrics(fitInfoML,simStack);
    validationInfo(n).fitInfoML = fitInfoML;
    
    % add info about individual parameters
    validationInfo(n).trueSpot1Fluo = trueInfo.GaussIntegral1;
    validationInfo(n).trueSpot2Fluo = trueInfo.GaussIntegral2;
    validationInfo(n).trueTotalFluo = trueInfo.GaussIntegralTot;

    validationInfo(n).bootSpot1Fluo = fitInfoBoot.GaussIntegral1;
    validationInfo(n).bootSpot2Fluo = fitInfoBoot.GaussIntegral2;
    validationInfo(n).bootTotalFluo = fitInfoBoot.GaussIntegralTot;
    
    validationInfo(n).mlSpot1Fluo = fitInfoML.GaussIntegral1;
    validationInfo(n).mlSpot2Fluo = fitInfoML.GaussIntegral2;
    validationInfo(n).mlTotalFluo = fitInfoML.GaussIntegralTot;
            
    validationInfo(n).fitSpot1Fluo = fitInfoRaw.GaussIntegral1;
    validationInfo(n).fitSpot2Fluo = fitInfoRaw.GaussIntegral2;
    validationInfo(n).fitTotalFluo = fitInfoRaw.GaussIntegralTot;
  
    
    validationInfo(n).trueTotalFluoRaw = trueInfo.GaussIntegralRaw;

    validationInfo(n).bootTotalFluoRaw = fitInfoBoot.GaussIntegralRaw;
    
    validationInfo(n).mlTotalFluoRaw = fitInfoML.GaussIntegralRaw;
    
    validationInfo(n).fitTotalFluoRaw = fitInfoRaw.GaussIntegralRaw;
    

    validationInfo(n).trueSpot1Pos = trueInfo.Spot1Pos;
    validationInfo(n).trueSpot2Pos = trueInfo.Spot2Pos;
    validationInfo(n).trueCenterPos = trueInfo.SpotCentroid;
   
    validationInfo(n).bootSpot1Pos = fitInfoBoot.Spot1Pos;
    validationInfo(n).bootSpot2Pos = fitInfoBoot.Spot2Pos;
    validationInfo(n).bootCenterPos = fitInfoBoot.SpotCentroid;
    
    validationInfo(n).mlSpot1Pos = fitInfoML.Spot1Pos;
    validationInfo(n).mlSpot2Pos = fitInfoML.Spot2Pos;
    validationInfo(n).mlCenterPos = fitInfoML.SpotCentroid;
    
    validationInfo(n).fitSpot1Pos = fitInfoRaw.Spot1Pos;
    validationInfo(n).fitSpot2Pos = fitInfoRaw.Spot2Pos;
    validationInfo(n).fitCenterPos = fitInfoRaw.SpotCentroid;
end
    
disp('Done.')  
%% %%%%%%%%%%%%%%%%%%%%%%%%%% Examine fit results %%%%%%%%%%%%%%%%%%%%%%%%% 
close all
labelCell = {'x','y','z'};
dimFactor = [PixelsSizeXY PixelsSizeXY PixelsSizeZ];

truePosArray = vertcat(validationInfo.trueCenterPos);
fitPosArray = vertcat(validationInfo.fitCenterPos);
bootPosArray = vertcat(validationInfo.bootCenterPos);
mlPosArray = vertcat(validationInfo.mlCenterPos);

cmap1 = brewermap([],'Set2');

for i = 1:3
    centroidFig = figure;
    hold on
    
    scatter(truePosArray(:,i),fitPosArray(:,i),'MArkerFaceColor',cmap1(1,:),'MArkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.4)
    scatter(truePosArray(:,i),bootPosArray(:,i),'MArkerFaceColor',cmap1(2,:),'MArkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.4)
    scatter(truePosArray(:,i),mlPosArray(:,i),'MArkerFaceColor',cmap1(3,:),'MArkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.4)

    xlabel(['true ' labelCell{i} ' position (\mu m)'])
    ylabel(['inferred ' labelCell{i} ' position (\mu m)'])
    grid on
    legend()
    set(gca,'FontSize',14)
    legend('single fit','bootstrap average','best fit','Location','southeast')
end

%%
fluoIntFig = figure;
hold on
cmap1 = brewermap([],'Set2');
scatter([validationInfo.trueTotalFluo],[validationInfo.fitTotalFluo],'MArkerFaceColor',cmap1(1,:),'MArkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.2)
scatter([validationInfo.trueTotalFluo],[validationInfo.bootTotalFluo],'MArkerFaceColor',cmap1(2,:),'MArkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.2)
scatter([validationInfo.trueTotalFluo],[validationInfo.mlTotalFluo],'MArkerFaceColor',cmap1(3,:),'MArkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.2)

xlabel('true fluorescence (au)')
ylabel('inferred fluorescence (au)')
grid on
legend('single fit','bootstrap average','best fit','Location','southeast')
set(gca,'FontSize',14)


fluoRawFig = figure;
hold on
cmap1 = brewermap([],'Set2');
scatter([validationInfo.trueTotalFluoRaw],[validationInfo.fitTotalFluoRaw],'MArkerFaceColor',cmap1(1,:),'MArkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.2)
scatter([validationInfo.trueTotalFluoRaw],[validationInfo.bootTotalFluoRaw],'MArkerFaceColor',cmap1(2,:),'MArkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.2)
scatter([validationInfo.trueTotalFluoRaw],[validationInfo.mlTotalFluoRaw],'MArkerFaceColor',cmap1(3,:),'MArkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.2)

xlabel('true fluorescence (au)')
ylabel('inferred fluorescence (au)')
grid on
legend('single fit','bootstrap average','best fit','Location','southeast')
set(gca,'FontSize',14)



%%
close all
fluo_err = [validationInfo.trueTotalFluoRaw];%([validationInfo.trueTotalFluoRaw]-[validationInfo.fitTotalFluoRaw]).^2;
x_err_boot = abs(truePosArray(:,2)-bootPosArray(:,2));
x_err_ml = abs(truePosArray(:,2)-mlPosArray(:,2));
x_err = abs(truePosArray(:,2)-fitPosArray(:,2));

f_err_ml = abs([validationInfo.trueTotalFluoRaw]-[validationInfo.mlTotalFluoRaw]);
f_err = abs([validationInfo.trueTotalFluoRaw]-[validationInfo.fitTotalFluoRaw]);

figure;
hold on
% plot(linspace(0,1,nSpotsToSimulate),sort(x_err_boot))
plot(linspace(0,1,nSpotsToSimulate),sort(x_err_ml))
plot(linspace(0,1,nSpotsToSimulate),sort(x_err))
grid on

figure;
hold on
% plot(linspace(0,1,nSpotsToSimulate),sort(x_err_boot))
plot(linspace(0,1,nSpotsToSimulate),sort(f_err_ml))
plot(linspace(0,1,nSpotsToSimulate),sort(f_err))
grid on