clear
close all

% set some basic simulation parameters
snipSize = [29 29 12];
stackSize = [50 50 14];
nSpotsToSimulate = 1e3;
nBootstraps =2;

% load Spots data for reference set
Prefix = '2020-03-15-2xDl-Ven_hbP2P-mCh_Leica_Zoom2_7uW14uW_05';
liveExperiment = LiveExperiment(Prefix);
FrameInfo = getFrameInfo(liveExperiment);
Spots = getSpots(liveExperiment);

% generate path to save figures
MS2CodePath = liveExperiment.MS2CodePath;
slashes = strfind(MS2CodePath,'\');
LivemRNAPath = MS2CodePath(1:slashes(end));
FigurePath = [LivemRNAPath 'SpotFittingValidation' filesep];
mkdir(FigurePath);

% get pixel dimensions
PixelSizeXY = FrameInfo(1).PixelSize;
PixelSizeZ = FrameInfo(1).ZStep;

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
makeOffsetSnip2Spot = @(params) normrnd(0,.5*params(11),size(mesh_z)) + ...
                                params(11) + ...
                                params(12)*mesh_y + params(13)*mesh_x + params(14)*mesh_z + ...
                                params(15)*mesh_y.^2 + params(16)*mesh_x.^2 + params(17)*mesh_z.^2;        
    
spot3DFun2Spot = @(params) simulate3DGaussSymmetric(snipSize, params(spot1ParamIndices))...
                                    + simulate3DGaussSymmetric(snipSize, params(spot2ParamIndices)) ...
                                    + makeOffsetSnip2Spot(params);
  
%initialize structure to store validation results
validationInfo = struct;
validationInfo(nSpotsToSimulate) = struct;

% store paramter values

startParallelPool(24, 0, 1);
if true %size(paramDelta,2) == 14
    xyMin = paramVecLb(11)/snipSize(1).^2;
    zMin = paramVecLb(11)/snipSize(3).^2;
    paramDelta(15:16) = 2*xyMin;
    paramDelta(17) = 2*zMin;
    paramVecLb(:,15:16) = -xyMin;
    paramVecLb(:,17) = -zMin;
end

disp('Simulating spots and fitting...')
parfor n = 1:nSpotsToSimulate
    % draw parameter values
    trueParams = rand(size(paramDelta)).*paramDelta + paramVecLb;
    % ensure that spot 1 is "below" spot 2
    spot1Score = sum(trueParams(4:6).^2);
    spot2Score = sum(trueParams(8:10).^2);
    orderIndices = 1:length(trueParams);
    if spot1Score > spot2Score        
        orderIndices(spot1ParamIndices) = spot2ParamIndices;
        orderIndices(spot2ParamIndices) = spot1ParamIndices;
    end   
    trueParams = trueParams(orderIndices);
    
    validationInfo(n).trueParams = trueParams;
    
    % create simulated snip
    simStack = spot3DFun2Spot(trueParams);    
    validationInfo(n).simStack = simStack;
    
    % test out "vanilla" initialization
    fitInfoRaw = fit3DGaussians(validationInfo(n).simStack,liveExperiment.pixelSize_nm,liveExperiment.zStep_um*1000,[],2,[]);    
    
    trueInfo = fit3DGaussians(validationInfo(n).simStack,liveExperiment.pixelSize_nm,liveExperiment.zStep_um*1000,[],2,trueParams);

    % calculate average for results            
    validationInfo(n).fitInfoRaw = fitInfoRaw;
    
    % add info about individual parameters
    validationInfo(n).trueSpot1Fluo = trueInfo.GaussIntegral1;
    validationInfo(n).trueSpot2Fluo = trueInfo.GaussIntegral2;
    validationInfo(n).trueTotalFluo = trueInfo.GaussIntegralTot;

    validationInfo(n).fitSpot1Fluo = fitInfoRaw.GaussIntegral1;
    validationInfo(n).fitSpot2Fluo = fitInfoRaw.GaussIntegral2;
    validationInfo(n).fitTotalFluo = fitInfoRaw.GaussIntegralTot;
  
    
    validationInfo(n).trueTotalFluoRaw = trueInfo.GaussIntegralRaw;   

    validationInfo(n).fitTotalFluoRaw = fitInfoRaw.GaussIntegralRaw;
    

    validationInfo(n).trueSpot1Pos = trueInfo.Spot1Pos;
    validationInfo(n).trueSpot2Pos = trueInfo.Spot2Pos;
    validationInfo(n).trueCenterPos = trueInfo.SpotCentroid;   

    validationInfo(n).fitSpot1Pos = fitInfoRaw.Spot1Pos;
    validationInfo(n).fitSpot2Pos = fitInfoRaw.Spot2Pos;
    validationInfo(n).fitCenterPos = fitInfoRaw.SpotCentroid;
end
    
disp('Done.')  
%% %%%%%%%%%%%%%%%%%%%%%%%%%% Examine fit results %%%%%%%%%%%%%%%%%%%%%%%%% 
close all
labelCell = {'x','y','z'};
dimFactor = [PixelSizeXY PixelSizeXY PixelSizeZ];

truePosArray = vertcat(validationInfo.trueCenterPos)-snipSize/2;
fitPosArray = vertcat(validationInfo.fitCenterPos)-snipSize/2;
% mlPosArray = vertcat(validationInfo.mlCenterPos);

cmap1 = brewermap([],'Set2');


centroidFig = figure;
hold on

scatter(truePosArray(:,1),fitPosArray(:,1),'MArkerFaceColor',cmap1(1,:),'MArkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.4)    
scatter(truePosArray(:,2),fitPosArray(:,2),'MArkerFaceColor',cmap1(2,:),'MArkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.4)    
scatter(truePosArray(:,3),fitPosArray(:,3),'MArkerFaceColor',cmap1(3,:),'MArkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.4)    

xlabel(['true position (pixels)'])
ylabel(['inferred position (pixels)'])
grid on    
set(gca,'FontSize',14)
legend('y','x','z','Location','southeast')

saveas(centroidFig,[FigurePath 'position_validation.png'])


%%%%%%%%%%%%%%%%%% fluorescence scatter %%%%%%%%%%%%%%%%%%
fluoIntFig = figure;
hold on
cmap1 = brewermap([],'Set2');
scatter([validationInfo.trueTotalFluo],[validationInfo.fitTotalFluo],'MArkerFaceColor',cmap1(2,:),'MArkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.5)
scatter([validationInfo.trueTotalFluoRaw],[validationInfo.fitTotalFluoRaw],'MArkerFaceColor',cmap1(3,:),'MArkerEdgeColor','k','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.5)

xlabel('true fluorescence (au)')
ylabel('inferred fluorescence (au)')
grid on
legend('Gaussian integral','raw intensity','Location','southeast')
set(gca,'FontSize',14)
saveas(fluoIntFig,[FigurePath 'spot_integral_validation.png'])

%%%%%%%%%%%%%%%%%% plot rank-ordered errors %%%%%%%%%%%%%%%%%%
close all
fluo_err = [validationInfo.trueTotalFluoRaw];%([validationInfo.trueTotalFluoRaw]-[validationInfo.fitTotalFluoRaw]).^2;
pos_err = sqrt(sum((truePosArray-fitPosArray).^2,2));
f_err = abs([validationInfo.trueTotalFluo]-[validationInfo.fitTotalFluo]);
f_err_raw = abs([validationInfo.trueTotalFluoRaw]-[validationInfo.fitTotalFluoRaw]);

posErrFig = figure;
hold on
plot(linspace(0,100,nSpotsToSimulate),sort(pos_err))
grid on
xlabel('percentile')
ylabel('error in centroid position (pixels)')
set(gca,'FontSize',14)
saveas(posErrFig,[FigurePath 'position_error_plot.png'])

fluoErrFig = figure;
hold on
% plot(linspace(0,1,nSpotsToSimulate),sort(x_err_boot))
plot(linspace(0,100,nSpotsToSimulate),sort(f_err))
plot(linspace(0,100,nSpotsToSimulate),sort(f_err_raw))
xlabel('percentile')
ylabel('error in raw spot intensity (au)')
set(gca,'FontSize',14)
grid on
legend('integral','raw intensity')
saveas(fluoErrFig,[FigurePath 'fluo_error_plot.png'])

