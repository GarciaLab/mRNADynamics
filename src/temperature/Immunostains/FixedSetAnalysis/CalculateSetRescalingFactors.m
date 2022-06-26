function AllCompiledEmbryos = CalculateSetRescalingFactors(AllCompiledEmbryos)
%%

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';
MasterSetPath = 'S:/Gabriella/Dropbox/ProteinProfiles/25CMasterSets.mat';
load(MasterSetPath, 'CombinedMean', 'CombinedSE', 'CombinedCounts', 'Slopes', 'Intercepts', 'Fits', 'SubsetsIncluded',...
    'Ts', 'Reps', 'CTstrings', 'SubsetsIncluded');
AllSetInfo = GetFixedSetPrefixInfo;
 
NumSets = length(AllSetInfo.Temperatures);

NChannels = 5;


RefSet = 14;
APbins = 0:0.025:1;


ControlFittingBins = zeros(NChannels, length(APbins), 'logical');
ControlFittingBins(2,:) = APbins >= 0.2 & APbins <= 0.8; % Hoechst
ControlFittingBins(3,:) = APbins >= 0.1 & APbins <= 0.9 ;%APbins >= 0.1 & APbins <= 0.5 ; % Bicoid
ControlFittingBins(4,:) = APbins >= 0.5 & APbins <= 0.9 ;  %APbins >= 0.5 & APbins <= 0.9 ;  % Knirps
ControlFittingBins(5,:) = APbins >= 0.1 & APbins <= 0.9 ;% (APbins >= 0.1 & APbins <= 0.5) | (APbins >= 0.7 & APbins <= 0.9); % Hunchback

ControlFittingBins = zeros(NChannels, NumAPbins, 'logical');
ControlFittingBins(2,:) = APbins >= 0.2 &  APbins <= 0.8 ;
ControlFittingBins(3,:) = APbins >= 0.1 &  APbins <= 0.5 ;
ControlFittingBins(4,:) = APbins >= 0.5 &  APbins <= 0.9 ;
ControlFittingBins(5,:) = (APbins >= 0.2 &  APbins <= 0.5) | (APbins >= 0.7 * APbins <= 0.9) ;
%%
if ~exist('AllCompiledEmbryos', 'var')
AllCompiledEmbryos = {};

for exp_index = 1:length(AllSetInfo.Temperatures)
if AllSetInfo.Flipped(exp_index)
    FlipString = 'yes';
else
    FlipString = 'no';
end
disp(['T = ', num2str(AllSetInfo.Temperatures(exp_index)), ', Rep: ', num2str(AllSetInfo.Replicates(exp_index)), ', Flipped: ', FlipString])
SetLabel = AllSetInfo.SetLabels{exp_index};
SetPrefixes = AllSetInfo.Prefixes{exp_index};
ProfFigPath = [AllSetsProfFigPath, SetLabel];
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
if ~isdir(ProfFigPath)
    mkdir(ProfFigPath)
end
if ~isdir(OutEmbryoPath)
    mkdir(OutEmbryoPath)
end
liveExperiments = cell(1, length(SetPrefixes));
for i = 1:length(SetPrefixes)
    liveExperiments{i} = LiveExperiment(SetPrefixes{i});
end
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
load(CEoutpath, 'CompiledEmbryos');
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
end
end



%%
warning('off','stats:LinearModel:RankDefDesignMat')
exp_idx = 1:NumSets;
AllSets = 1:NumSets;
counts_1sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_1sigma, 2));
counts_2sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_2sigma, 2));
ControlledProfiles = NaN([size(AllCompiledEmbryos{1}.FixCorrectedSmoothedAvgAPProfiles.Control), NumSets]);
for i = 1:NumSets
counts_1sigma(i,:) = AllCompiledEmbryos{exp_idx(i)}.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_1sigma;
counts_2sigma(i,:) = AllCompiledEmbryos{exp_idx(i)}.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_2sigma;
ControlledProfiles(:,:,:,i) = AllCompiledEmbryos{exp_idx(i)}.FixCorrectedSmoothedAvgAPProfiles.Control;
end
FitInfo = cell(NumSets, NumSets, NChannels);
ScalingFactors = NaN(NumSets, NumSets, NChannels);
ScalingSEs = NaN(NumSets, NumSets, NChannels);
ScalingIntercepts = NaN(NumSets, NumSets, NChannels);
ScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
GoodOverlapIndices = cell(NumSets, NumSets);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
OkScalingIntercepts = NaN(NumSets, NumSets, NChannels);
OkScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
OkOverlapIndices = cell(NumSets, NumSets);
for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        OkOverlapIndices{exp2, exp1} = (counts_1sigma(i,:)+counts_2sigma(i,:) >= 5) & (counts_1sigma(j,:)+counts_2sigma(j,:) >= 5);
        GoodOverlapIndices{exp2, exp1} = (counts_1sigma(i,:) >= 5) & (counts_1sigma(j,:) >= 5);
        if sum(GoodOverlapIndices{exp2, exp1}) >= 8
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( GoodOverlapIndices{exp2, exp1},:);
                FitSet1 = FitSet1(:);
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
                FitSet2 = FitSet2( GoodOverlapIndices{exp2, exp1},:);
                FitSet2 = FitSet2(:);
                dlm = fitlm(FitSet2, FitSet1);
                FitInfo{exp2, exp1, ch_index} = dlm;
                ScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(2);
                ScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(2);
                ScalingIntercepts(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(1);
                ScalingInterceptSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(1);
            end
        end
        
        if sum(OkOverlapIndices{exp2, exp1}) >= 8
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1},:);
                FitSet1 = FitSet1(:);
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
                FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1},:);
                FitSet2 = FitSet2(:);
                dlm = fitlm(FitSet2, FitSet1);
                OkFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(2);
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(2);
                OkScalingIntercepts(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(1);
                OkScalingInterceptSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(1);
            end
        end
    end
end

for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            ScalingFactors(exp1, exp2, ch_index) = 1/ScalingFactors(exp2, exp1, ch_index);
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/(ScalingFactors(exp2, exp1, ch_index)^2);
            ScalingIntercepts(exp1, exp2, ch_index) = -ScalingIntercepts(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            ScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(ScalingInterceptSEs(exp2, exp1, ch_index)^2+ScalingIntercepts(exp1, exp2, ch_index)^2*ScalingSEs(exp2, exp1, ch_index)^2)/ScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingIntercepts(exp1, exp2, ch_index) = -OkScalingIntercepts(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(OkScalingInterceptSEs(exp2, exp1, ch_index)^2+OkScalingIntercepts(exp1, exp2, ch_index)^2*OkScalingSEs(exp2, exp1, ch_index)^2)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingIntercepts = NaN(NumSets, NChannels);
SetScalingInterceptSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSet, 2:NChannels) = 1;
SetScalingSEs(RefSet, 2:NChannels) = 0;
SetScalingIntercepts(RefSet, 2:NChannels) = 0;
SetScalingInterceptSEs(RefSet, 2:NChannels) = 0;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, RefSet))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSet, ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSet, ch_index);
        SetScalingIntercepts(i, ch_index) = ScalingIntercepts(i, RefSet, ch_index);
        SetScalingInterceptSEs(i, ch_index) = ScalingInterceptSEs(i, RefSet, ch_index);
        if isnan(SetScalingFactors(i, ch_index) )
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSet, ch_index);
            SetScalingIntercepts(i, ch_index) = OkScalingIntercepts(i, RefSet, ch_index);
            SetScalingInterceptSEs(i, ch_index) = OkScalingInterceptSEs(i, RefSet, ch_index);
        end

    end
end

for exp_index = 1:length(AllSetInfo.Temperatures)
SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
CompiledEmbryos = AllCompiledEmbryos{exp_index};
CompiledEmbryos.FixCorrectedControlScalingFactors = SetScalingFactors(exp_index,:);
CompiledEmbryos.FixCorrectedControlScalingSEs = SetScalingSEs(exp_index,:);
CompiledEmbryos.FixCorrectedControlScalingIntercepts = SetScalingIntercepts(exp_index,:);
CompiledEmbryos.FixCorrectedControlScalingInterceptSEs = SetScalingInterceptSEs(exp_index,:);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
save(CEoutpath, 'CompiledEmbryos');
end

%% Fit Universal Scaling using Dubuis Embryo Times
exp_idx = 1:NumSets;
AllSets = 1:NumSets;
counts_1sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_1sigma, 2));
counts_2sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_2sigma, 2));
counts_3sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_2sigma, 2));
ControlledProfiles = NaN([size(AllCompiledEmbryos{1}.DubuisTimeSmoothedAvgAPProfiles.Control), NumSets]);
for i = 1:NumSets
counts_1sigma(i,:) = AllCompiledEmbryos{exp_idx(i)}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_1sigma;
counts_2sigma(i,:) = AllCompiledEmbryos{exp_idx(i)}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_2sigma;
counts_3sigma(i,:) = AllCompiledEmbryos{exp_idx(i)}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_3sigma;
ControlledProfiles(:,:,:,i) = AllCompiledEmbryos{exp_idx(i)}.DubuisTimeSmoothedAvgAPProfiles.Control;
end
FitInfo = cell(NumSets, NumSets, NChannels);
ScalingFactors = NaN(NumSets, NumSets, NChannels);
ScalingSEs = NaN(NumSets, NumSets, NChannels);
ScalingIntercepts = NaN(NumSets, NumSets, NChannels);
ScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
GoodOverlapIndices = cell(NumSets, NumSets);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
OkScalingIntercepts = NaN(NumSets, NumSets, NChannels);
OkScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
OkOverlapIndices = cell(NumSets, NumSets);
for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        OkOverlapIndices{exp2, exp1} = (counts_1sigma(i,:)+counts_2sigma(i,:) >= 5) & (counts_1sigma(j,:)+counts_2sigma(j,:) >= 5);
        GoodOverlapIndices{exp2, exp1} = (counts_1sigma(i,:) >= 5) & (counts_1sigma(j,:) >= 5);
        if sum(GoodOverlapIndices{exp2, exp1}) >= 8
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( GoodOverlapIndices{exp2, exp1},:);
                FitSet1 = FitSet1(:);
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
                FitSet2 = FitSet2( GoodOverlapIndices{exp2, exp1},:);
                FitSet2 = FitSet2(:);
                dlm = fitlm(FitSet2, FitSet1);
                FitInfo{exp2, exp1, ch_index} = dlm;
                ScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(2);
                ScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(2);
                ScalingIntercepts(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(1);
                ScalingInterceptSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(1);
            end

        end
        
        if sum(OkOverlapIndices{exp2, exp1}) >= 8
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1},:);
                FitSet1 = FitSet1(:);
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
                FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1},:);
                FitSet2 = FitSet2(:);
                dlm = fitlm(FitSet2, FitSet1);
                OkFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(2);
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(2);
                OkScalingIntercepts(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(1);
                OkScalingInterceptSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(1);
            end
        end
    end
end

for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            ScalingFactors(exp1, exp2, ch_index) = 1/ScalingFactors(exp2, exp1, ch_index);
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/(ScalingFactors(exp2, exp1, ch_index)^2);
            ScalingIntercepts(exp1, exp2, ch_index) = -ScalingIntercepts(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            ScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(ScalingInterceptSEs(exp2, exp1, ch_index)^2+ScalingIntercepts(exp1, exp2, ch_index)^2*ScalingSEs(exp2, exp1, ch_index)^2)/ScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingIntercepts(exp1, exp2, ch_index) = -OkScalingIntercepts(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(OkScalingInterceptSEs(exp2, exp1, ch_index)^2+OkScalingIntercepts(exp1, exp2, ch_index)^2*OkScalingSEs(exp2, exp1, ch_index)^2)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingIntercepts = NaN(NumSets, NChannels);
SetScalingInterceptSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSet, 2:NChannels) = 1;
SetScalingSEs(RefSet, 2:NChannels) = 0;
SetScalingIntercepts(RefSet, 2:NChannels) = 0;
SetScalingInterceptSEs(RefSet, 2:NChannels) = 0;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, RefSet))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSet, ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSet, ch_index);
        SetScalingIntercepts(i, ch_index) = ScalingIntercepts(i, RefSet, ch_index);
        SetScalingInterceptSEs(i, ch_index) = ScalingInterceptSEs(i, RefSet, ch_index);
        if isnan(SetScalingFactors(i, ch_index) )
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSet, ch_index);
            SetScalingIntercepts(i, ch_index) = OkScalingIntercepts(i, RefSet, ch_index);
            SetScalingInterceptSEs(i, ch_index) = OkScalingInterceptSEs(i, RefSet, ch_index);
        end

    end
end

for exp_index = 1:length(AllSetInfo.Temperatures)
SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
CompiledEmbryos = AllCompiledEmbryos{exp_index};
CompiledEmbryos.DubuisTimesControlScalingFactors = SetScalingFactors(exp_index,:);
CompiledEmbryos.DubuisTimesControlScalingSEs = SetScalingSEs(exp_index,:);
CompiledEmbryos.DubuisTimesControlScalingIntercepts = SetScalingIntercepts(exp_index,:);
CompiledEmbryos.DubuisTimesControlScalingInterceptSEs = SetScalingInterceptSEs(exp_index,:);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
save(CEoutpath, 'CompiledEmbryos');
end



%% Fit Universal Scaling using HisRFP 25C Times
exp_idx = 1:NumSets;
AllSets = 1:NumSets;
counts_1sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.HisRFPTimesSmoothedProfiles.ControlCounts.counts_within_1sigma, 2));
counts_2sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.HisRFPTimesSmoothedProfiles.ControlCounts.counts_within_2sigma, 2));
counts_3sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.HisRFPTimesSmoothedProfiles.ControlCounts.counts_within_2sigma, 2));
ControlledProfiles = NaN([size(AllCompiledEmbryos{1}.HisRFP25CTimeSmoothedAvgAPProfiles.Control), NumSets]);
for i = 1:NumSets
counts_1sigma(i,:) = AllCompiledEmbryos{exp_idx(i)}.HisRFPTimesSmoothedProfiles.ControlCounts.counts_within_1sigma;
counts_2sigma(i,:) = AllCompiledEmbryos{exp_idx(i)}.HisRFPTimesSmoothedProfiles.ControlCounts.counts_within_2sigma;
counts_3sigma(i,:) = AllCompiledEmbryos{exp_idx(i)}.HisRFPTimesSmoothedProfiles.ControlCounts.counts_within_3sigma;
ControlledProfiles(:,:,:,i) = AllCompiledEmbryos{exp_idx(i)}.HisRFP25CTimeSmoothedAvgAPProfiles.Control;
end
FitInfo = cell(NumSets, NumSets, NChannels);
ScalingFactors = NaN(NumSets, NumSets, NChannels);
ScalingSEs = NaN(NumSets, NumSets, NChannels);
ScalingIntercepts = NaN(NumSets, NumSets, NChannels);
ScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
GoodOverlapIndices = cell(NumSets, NumSets);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
OkScalingIntercepts = NaN(NumSets, NumSets, NChannels);
OkScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
OkOverlapIndices = cell(NumSets, NumSets);


for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        OkOverlapIndices{exp2, exp1} = (counts_1sigma(i,:)+counts_2sigma(i,:)+counts_3sigma(i,:) >= 5) & (counts_1sigma(j,:)+counts_2sigma(j,:)+counts_3sigma(j,:) >= 5);
        GoodOverlapIndices{exp2, exp1} = (counts_1sigma(i,:)+counts_2sigma(i,:) >= 5) & (counts_1sigma(j,:)+counts_2sigma(j,:) >= 5);
        if sum(GoodOverlapIndices{exp2, exp1}) >= 8
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( GoodOverlapIndices{exp2, exp1},:);
                FitSet1 = FitSet1(:);
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
                FitSet2 = FitSet2( GoodOverlapIndices{exp2, exp1},:);
                FitSet2 = FitSet2(:);
                dlm = fitlm(FitSet2, FitSet1);
                FitInfo{exp2, exp1, ch_index} = dlm;
                ScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(2);
                ScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(2);
                ScalingIntercepts(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(1);
                ScalingInterceptSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(1);
            end
        end
        
        if sum(OkOverlapIndices{exp2, exp1}) >= 8
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1},:);
                FitSet1 = FitSet1(:);
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
                FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1},:);
                FitSet2 = FitSet2(:);
                dlm = fitlm(FitSet2, FitSet1);
                OkFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(2);
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(2);
                OkScalingIntercepts(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(1);
                OkScalingInterceptSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(1);
            end
        end
    end
end

for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            ScalingFactors(exp1, exp2, ch_index) = 1/ScalingFactors(exp2, exp1, ch_index);
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/(ScalingFactors(exp2, exp1, ch_index)^2);
            ScalingIntercepts(exp1, exp2, ch_index) = -ScalingIntercepts(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            ScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(ScalingInterceptSEs(exp2, exp1, ch_index)^2+ScalingIntercepts(exp1, exp2, ch_index)^2*ScalingSEs(exp2, exp1, ch_index)^2)/ScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingIntercepts(exp1, exp2, ch_index) = -OkScalingIntercepts(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(OkScalingInterceptSEs(exp2, exp1, ch_index)^2+OkScalingIntercepts(exp1, exp2, ch_index)^2*OkScalingSEs(exp2, exp1, ch_index)^2)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingIntercepts = NaN(NumSets, NChannels);
SetScalingInterceptSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSet, 2:NChannels) = 1;
SetScalingSEs(RefSet, 2:NChannels) = 0;
SetScalingIntercepts(RefSet, 2:NChannels) = 0;
SetScalingInterceptSEs(RefSet, 2:NChannels) = 0;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, RefSet))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSet, ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSet, ch_index);
        SetScalingIntercepts(i, ch_index) = ScalingIntercepts(i, RefSet, ch_index);
        SetScalingInterceptSEs(i, ch_index) = ScalingInterceptSEs(i, RefSet, ch_index);
        if isnan(SetScalingFactors(i, ch_index) )
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSet, ch_index);
            SetScalingIntercepts(i, ch_index) = OkScalingIntercepts(i, RefSet, ch_index);
            SetScalingInterceptSEs(i, ch_index) = OkScalingInterceptSEs(i, RefSet, ch_index);
        end

    end
end

for exp_index = 1:length(AllSetInfo.Temperatures)
SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
CompiledEmbryos = AllCompiledEmbryos{exp_index};
CompiledEmbryos.HisRFP25CTimesControlScalingFactors = SetScalingFactors(exp_index,:);
CompiledEmbryos.HisRFP25CTimesControlScalingSEs = SetScalingSEs(exp_index,:);
CompiledEmbryos.HisRFP25CTimesControlScalingIntercepts = SetScalingIntercepts(exp_index,:);
CompiledEmbryos.HisRFP25CTimesControlScalingInterceptSEs = SetScalingInterceptSEs(exp_index,:);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
save(CEoutpath, 'CompiledEmbryos');
end

%% Windowed Dubuis Time Profile Info 

warning('off','stats:LinearModel:RankDefDesignMat')
exp_idx = 1:NumSets;
AllSets = 1:NumSets;
CountDims = [size(AllCompiledEmbryos{1}.WindowedProfiles.DubuisTime.Control.count, 1),...
    size(AllCompiledEmbryos{1}.WindowedProfiles.DubuisTime.Control.count, 2), NumSets];
counts = zeros(CountDims);

ControlledProfiles = NaN([size(AllCompiledEmbryos{1}.WindowedProfiles.DubuisTime.Control.count) NumSets]);
for i = 1:NumSets
counts(:,:,i) = AllCompiledEmbryos{exp_idx(i)}.WindowedProfiles.DubuisTime.Control.count(:,:,3);
ControlledProfiles(:,:,:,i) = AllCompiledEmbryos{exp_idx(i)}.WindowedProfiles.DubuisTime.Control.mean;
end
FitInfo = cell(NumSets, NumSets, NChannels);
ScalingFactors = NaN(NumSets, NumSets, NChannels);
ScalingSEs = NaN(NumSets, NumSets, NChannels);
ScalingIntercepts = NaN(NumSets, NumSets, NChannels);
ScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
BestOverlapIndices = cell(NumSets, NumSets, NChannels);

GoodFitInfo = cell(NumSets, NumSets, NChannels);
GoodScalingFactors = NaN(NumSets, NumSets, NChannels);
GoodScalingSEs = NaN(NumSets, NumSets, NChannels);
GoodScalingIntercepts = NaN(NumSets, NumSets, NChannels);
GoodScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
GoodOverlapIndices = cell(NumSets, NumSets, NChannels);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
OkScalingIntercepts = NaN(NumSets, NumSets, NChannels);
OkScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
OkOverlapIndices = cell(NumSets, NumSets, NChannels);
for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            BestOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 5) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 5);
            GoodOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 3) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 3);
            OkOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 2) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 2);
            if sum(sum((OkOverlapIndices{exp2, exp1, ch_index} )) >= 1) >= 1
                
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1, ch_index}).';
        
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
               FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1, ch_index}).';
       
                dlm = fitlm(FitSet2, FitSet1);%,'Weights',  1./(FitSD1.^2));
                OkFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(2);
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(2);
                OkScalingIntercepts(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(1);
                OkScalingInterceptSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(1);

            end
        end
        
    end
end


for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            ScalingFactors(exp1, exp2, ch_index) = 1/ScalingFactors(exp2, exp1, ch_index);
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/(ScalingFactors(exp2, exp1, ch_index)^2);
            ScalingIntercepts(exp1, exp2, ch_index) = -ScalingIntercepts(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            ScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(ScalingInterceptSEs(exp2, exp1, ch_index)^2+ScalingIntercepts(exp1, exp2, ch_index)^2*ScalingSEs(exp2, exp1, ch_index)^2)/ScalingFactors(exp2, exp1, ch_index);
            GoodScalingFactors(exp1, exp2, ch_index) = 1/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingSEs(exp1, exp2, ch_index) = GoodScalingSEs(exp2, exp1, ch_index)/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingIntercepts(exp1, exp2, ch_index) = -GoodScalingIntercepts(exp2, exp1, ch_index)/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(GoodScalingInterceptSEs(exp2, exp1, ch_index)^2+GoodScalingIntercepts(exp1, exp2, ch_index)^2*GoodScalingSEs(exp2, exp1, ch_index)^2)/GoodScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingIntercepts(exp1, exp2, ch_index) = -OkScalingIntercepts(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(OkScalingInterceptSEs(exp2, exp1, ch_index)^2+OkScalingIntercepts(exp1, exp2, ch_index)^2*OkScalingSEs(exp2, exp1, ch_index)^2)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingIntercepts = NaN(NumSets, NChannels);
SetScalingInterceptSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSet, 2:NChannels) = 1;
SetScalingSEs(RefSet, 2:NChannels) = 0;
SetScalingIntercepts(RefSet, 2:NChannels) = 0;
SetScalingInterceptSEs(RefSet, 2:NChannels) = 0;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, RefSet))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSet, ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSet, ch_index);
        SetScalingIntercepts(i, ch_index) = ScalingIntercepts(i, RefSet, ch_index);
        SetScalingInterceptSEs(i, ch_index) = ScalingInterceptSEs(i, RefSet, ch_index);
        if isnan(SetScalingFactors(i, ch_index) )
            SetScalingFactors(i, ch_index) = GoodScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = GoodScalingSEs(i, RefSet, ch_index);
            SetScalingIntercepts(i, ch_index) = GoodScalingIntercepts(i, RefSet, ch_index);
            SetScalingInterceptSEs(i, ch_index) = GoodScalingInterceptSEs(i, RefSet, ch_index);
        end
        if isnan(SetScalingFactors(i, ch_index) )
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSet, ch_index);
            SetScalingIntercepts(i, ch_index) = OkScalingIntercepts(i, RefSet, ch_index);
            SetScalingInterceptSEs(i, ch_index) = OkScalingInterceptSEs(i, RefSet, ch_index);
        end

    end
end

for exp_index = 1:length(AllSetInfo.Temperatures)
SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
CompiledEmbryos = AllCompiledEmbryos{exp_index};
CompiledEmbryos.WindowedDubuisTimesControlScalingFactors = SetScalingFactors(exp_index,:);
CompiledEmbryos.WindowedDubuisTimesControlScalingSEs = SetScalingSEs(exp_index,:);
CompiledEmbryos.WindowedDubuisTimesControlScalingIntercepts = SetScalingIntercepts(exp_index,:);
CompiledEmbryos.WindowedDubuisTimesControlScalingInterceptSEs = SetScalingInterceptSEs(exp_index,:);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
save(CEoutpath, 'CompiledEmbryos');
end

%% Windowed Delta FC Profile Info 

warning('off','stats:LinearModel:RankDefDesignMat')
RefSet = 14;
exp_idx = 1:NumSets;
AllSets = 1:NumSets;
CountDims = [size(AllCompiledEmbryos{1}.WindowedProfiles.DeltaFC.Control.count, 1),...
    size(AllCompiledEmbryos{1}.WindowedProfiles.DeltaFC.Control.count, 2), NumSets];
counts = zeros(CountDims);

ControlledProfiles = NaN([size(AllCompiledEmbryos{1}.WindowedProfiles.DeltaFC.Control.count) NumSets]);
ControlledProfileSDs = NaN([size(AllCompiledEmbryos{1}.WindowedProfiles.DeltaFC.Control.count) NumSets]);
for i = 1:NumSets
counts(:,:,i) = AllCompiledEmbryos{exp_idx(i)}.WindowedProfiles.DeltaFC.Control.count(:,:,3);
ControlledProfiles(:,:,:,i) = AllCompiledEmbryos{exp_idx(i)}.WindowedProfiles.DeltaFC.Control.mean;
ControlledProfileSDs(:,:,:,i) = AllCompiledEmbryos{exp_idx(i)}.WindowedProfiles.DeltaFC.Control.std;
end

FitInfo = cell(NumSets, NumSets, NChannels);
ScalingFactors = NaN(NumSets, NumSets, NChannels);
ScalingSEs = NaN(NumSets, NumSets, NChannels);
ScalingIntercepts = NaN(NumSets, NumSets, NChannels);
ScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
BestOverlapIndices = cell(NumSets, NumSets, NChannels);

GoodFitInfo = cell(NumSets, NumSets, NChannels);
GoodScalingFactors = NaN(NumSets, NumSets, NChannels);
GoodScalingSEs = NaN(NumSets, NumSets, NChannels);
GoodScalingIntercepts = NaN(NumSets, NumSets, NChannels);
GoodScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
GoodOverlapIndices = cell(NumSets, NumSets, NChannels);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
OkScalingIntercepts = NaN(NumSets, NumSets, NChannels);
OkScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
OkOverlapIndices = cell(NumSets, NumSets, NChannels);
for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            BestOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 5) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 5);
            GoodOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 3) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 3);
            OkOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 2) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 2);
%             if sum(sum((BestOverlapIndices{exp2, exp1, ch_index} )) >= 5) >= 1
%                 
%                 FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
%                 FitSE1 =  ControlledProfileSEs(:,ControlFittingBins(ch_index, :),ch_index,i);
%                 FitSet1 = FitSet1( BestOverlapIndices{exp2, exp1, ch_index}).';
%                 FitSE1 = FitSE1(    BestOverlapIndices{exp2, exp1, ch_index}).';
%                 FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
%                FitSet2 = FitSet2( BestOverlapIndices{exp2, exp1, ch_index}).';
%                 FitSE2 =  ControlledProfileSEs(:,ControlFittingBins(ch_index, :),ch_index,j);
%                 FitSE2 = FitSE2( BestOverlapIndices{exp2, exp1, ch_index}).';
%                
%                 dlm = fitlm(FitSet2, FitSet1);
%                 FitInfo{exp2, exp1, ch_index} = dlm;
%                 ScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(2);
%                 ScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(2);
%                 ScalingIntercepts(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(1);
%                 ScalingInterceptSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(1);
%                 
%             end
%             if sum(sum((GoodOverlapIndices{exp2, exp1, ch_index} )) >= 5) >= 1
%                 
%                 FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
%                 FitSet1 = FitSet1( GoodOverlapIndices{exp2, exp1, ch_index}).';
%                 FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
%                 FitSet2 = FitSet2( GoodOverlapIndices{exp2, exp1, ch_index}).';
%                
%                 dlm = fitlm(FitSet2, FitSet1);
%                 GoodFitInfo{exp2, exp1, ch_index} = dlm;
%                 GoodScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(2);
%                 GoodScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(2);
%                 GoodScalingIntercepts(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(1);
%                 GoodScalingInterceptSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(1);
%                 
            
            if sum(sum((OkOverlapIndices{exp2, exp1, ch_index} )) >= 1) >= 1
                
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1, ch_index}).';
        
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
               FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1, ch_index}).';
       
                dlm = fitlm(FitSet2, FitSet1);%,'Weights',  1./(FitSD1.^2));
                OkFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(2);
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(2);
                OkScalingIntercepts(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(1);
                OkScalingInterceptSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(1);

            end
        end
        
    end
end


for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            ScalingFactors(exp1, exp2, ch_index) = 1/ScalingFactors(exp2, exp1, ch_index);
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/(ScalingFactors(exp2, exp1, ch_index)^2);
            ScalingIntercepts(exp1, exp2, ch_index) = -ScalingIntercepts(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            ScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(ScalingInterceptSEs(exp2, exp1, ch_index)^2+ScalingIntercepts(exp1, exp2, ch_index)^2*ScalingSEs(exp2, exp1, ch_index)^2)/ScalingFactors(exp2, exp1, ch_index);
            GoodScalingFactors(exp1, exp2, ch_index) = 1/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingSEs(exp1, exp2, ch_index) = GoodScalingSEs(exp2, exp1, ch_index)/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingIntercepts(exp1, exp2, ch_index) = -GoodScalingIntercepts(exp2, exp1, ch_index)/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(GoodScalingInterceptSEs(exp2, exp1, ch_index)^2+GoodScalingIntercepts(exp1, exp2, ch_index)^2*GoodScalingSEs(exp2, exp1, ch_index)^2)/GoodScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingIntercepts(exp1, exp2, ch_index) = -OkScalingIntercepts(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(OkScalingInterceptSEs(exp2, exp1, ch_index)^2+OkScalingIntercepts(exp1, exp2, ch_index)^2*OkScalingSEs(exp2, exp1, ch_index)^2)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingIntercepts = NaN(NumSets, NChannels);
SetScalingInterceptSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSet, 2:NChannels) = 1;
SetScalingSEs(RefSet, 2:NChannels) = 0;
SetScalingIntercepts(RefSet, 2:NChannels) = 0;
SetScalingInterceptSEs(RefSet, 2:NChannels) = 0;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, RefSet))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSet, ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSet, ch_index);
        SetScalingIntercepts(i, ch_index) = ScalingIntercepts(i, RefSet, ch_index);
        SetScalingInterceptSEs(i, ch_index) = ScalingInterceptSEs(i, RefSet, ch_index);
        if isnan(SetScalingFactors(i, ch_index) )
            SetScalingFactors(i, ch_index) = GoodScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = GoodScalingSEs(i, RefSet, ch_index);
            SetScalingIntercepts(i, ch_index) = GoodScalingIntercepts(i, RefSet, ch_index);
            SetScalingInterceptSEs(i, ch_index) = GoodScalingInterceptSEs(i, RefSet, ch_index);
        end
        if isnan(SetScalingFactors(i, ch_index) )
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSet, ch_index);
            SetScalingIntercepts(i, ch_index) = OkScalingIntercepts(i, RefSet, ch_index);
            SetScalingInterceptSEs(i, ch_index) = OkScalingInterceptSEs(i, RefSet, ch_index);
        end

    end
end

for exp_index = 1:length(AllSetInfo.Temperatures)
SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
CompiledEmbryos = AllCompiledEmbryos{exp_index};
CompiledEmbryos.WindowedDeltaFCControlScalingFactors = SetScalingFactors(exp_index,:);
CompiledEmbryos.WindowedDeltaFCControlScalingSEs = SetScalingSEs(exp_index,:);
CompiledEmbryos.WindowedDeltaFCControlScalingIntercepts = SetScalingIntercepts(exp_index,:);
CompiledEmbryos.WindowedDeltaFCControlScalingInterceptSEs = SetScalingInterceptSEs(exp_index,:);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
save(CEoutpath, 'CompiledEmbryos');
end



%% Windowed yw25C Profile Info 

warning('off','stats:LinearModel:RankDefDesignMat')
RefSet = 14;
exp_idx = 1:NumSets;
AllSets = 1:NumSets;
CountDims = [size(AllCompiledEmbryos{1}.WindowedProfiles.yw25CTime.Control.count, 1),...
    size(AllCompiledEmbryos{1}.WindowedProfiles.yw25CTime.Control.count, 2), NumSets];
counts = zeros(CountDims);

ControlledProfiles = NaN([size(AllCompiledEmbryos{1}.WindowedProfiles.yw25CTime.Control.count) NumSets]);
for i = 1:NumSets
counts(:,:,i) = AllCompiledEmbryos{exp_idx(i)}.WindowedProfiles.yw25CTime.Control.count(:,:,3);
ControlledProfiles(:,:,:,i) = AllCompiledEmbryos{exp_idx(i)}.WindowedProfiles.yw25CTime.Control.mean;
end

FitInfo = cell(NumSets, NumSets, NChannels);
ScalingFactors = NaN(NumSets, NumSets, NChannels);
ScalingSEs = NaN(NumSets, NumSets, NChannels);
ScalingIntercepts = NaN(NumSets, NumSets, NChannels);
ScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
BestOverlapIndices = cell(NumSets, NumSets, NChannels);

GoodFitInfo = cell(NumSets, NumSets, NChannels);
GoodScalingFactors = NaN(NumSets, NumSets, NChannels);
GoodScalingSEs = NaN(NumSets, NumSets, NChannels);
GoodScalingIntercepts = NaN(NumSets, NumSets, NChannels);
GoodScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
GoodOverlapIndices = cell(NumSets, NumSets, NChannels);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
OkScalingIntercepts = NaN(NumSets, NumSets, NChannels);
OkScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
OkOverlapIndices = cell(NumSets, NumSets, NChannels);
for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            BestOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 5) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 5);
            GoodOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 3) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 3);
            OkOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 2) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 2);
            if sum(sum((OkOverlapIndices{exp2, exp1, ch_index} )) >= 1) >= 1
                
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1, ch_index}).';
        
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
               FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1, ch_index}).';
       
                dlm = fitlm(FitSet2, FitSet1);%,'Weights',  1./(FitSD1.^2));
                OkFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(2);
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(2);
                OkScalingIntercepts(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(1);
                OkScalingInterceptSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(1);

            end
        end
        
    end
end


for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            ScalingFactors(exp1, exp2, ch_index) = 1/ScalingFactors(exp2, exp1, ch_index);
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/(ScalingFactors(exp2, exp1, ch_index)^2);
            ScalingIntercepts(exp1, exp2, ch_index) = -ScalingIntercepts(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            ScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(ScalingInterceptSEs(exp2, exp1, ch_index)^2+ScalingIntercepts(exp1, exp2, ch_index)^2*ScalingSEs(exp2, exp1, ch_index)^2)/ScalingFactors(exp2, exp1, ch_index);
            GoodScalingFactors(exp1, exp2, ch_index) = 1/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingSEs(exp1, exp2, ch_index) = GoodScalingSEs(exp2, exp1, ch_index)/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingIntercepts(exp1, exp2, ch_index) = -GoodScalingIntercepts(exp2, exp1, ch_index)/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(GoodScalingInterceptSEs(exp2, exp1, ch_index)^2+GoodScalingIntercepts(exp1, exp2, ch_index)^2*GoodScalingSEs(exp2, exp1, ch_index)^2)/GoodScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingIntercepts(exp1, exp2, ch_index) = -OkScalingIntercepts(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(OkScalingInterceptSEs(exp2, exp1, ch_index)^2+OkScalingIntercepts(exp1, exp2, ch_index)^2*OkScalingSEs(exp2, exp1, ch_index)^2)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingIntercepts = NaN(NumSets, NChannels);
SetScalingInterceptSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSet, 2:NChannels) = 1;
SetScalingSEs(RefSet, 2:NChannels) = 0;
SetScalingIntercepts(RefSet, 2:NChannels) = 0;
SetScalingInterceptSEs(RefSet, 2:NChannels) = 0;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, RefSet))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSet, ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSet, ch_index);
        SetScalingIntercepts(i, ch_index) = ScalingIntercepts(i, RefSet, ch_index);
        SetScalingInterceptSEs(i, ch_index) = ScalingInterceptSEs(i, RefSet, ch_index);
        if isnan(SetScalingFactors(i, ch_index) )
            SetScalingFactors(i, ch_index) = GoodScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = GoodScalingSEs(i, RefSet, ch_index);
            SetScalingIntercepts(i, ch_index) = GoodScalingIntercepts(i, RefSet, ch_index);
            SetScalingInterceptSEs(i, ch_index) = GoodScalingInterceptSEs(i, RefSet, ch_index);
        end
        if isnan(SetScalingFactors(i, ch_index) )
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSet, ch_index);
            SetScalingIntercepts(i, ch_index) = OkScalingIntercepts(i, RefSet, ch_index);
            SetScalingInterceptSEs(i, ch_index) = OkScalingInterceptSEs(i, RefSet, ch_index);
        end

    end
end

for exp_index = 1:length(AllSetInfo.Temperatures)
    SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
CompiledEmbryos = AllCompiledEmbryos{exp_index};
CompiledEmbryos.Windowedyw25CTimesControlScalingFactors = SetScalingFactors(exp_index,:);
CompiledEmbryos.Windowedyw25CTimesControlScalingSEs = SetScalingSEs(exp_index,:);
CompiledEmbryos.Windowedyw25CTimesControlScalingIntercepts = SetScalingIntercepts(exp_index,:);
CompiledEmbryos.Windowedyw25CTimesControlScalingInterceptSEs = SetScalingInterceptSEs(exp_index,:);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
save(CEoutpath, 'CompiledEmbryos');
end


%% Windowed HisRFP 25C Profile Info 

warning('off','stats:LinearModel:RankDefDesignMat')
RefSet = 14;
exp_idx = 1:NumSets;
AllSets = 1:NumSets;
CountDims = [size(AllCompiledEmbryos{1}.WindowedProfiles.HisRFP25CTime.Control.count, 1),...
    size(AllCompiledEmbryos{1}.WindowedProfiles.HisRFP25CTime.Control.count, 2), NumSets];
counts = zeros(CountDims);

ControlledProfiles = NaN([size(AllCompiledEmbryos{1}.WindowedProfiles.HisRFP25CTime.Control.count) NumSets]);
for i = 1:NumSets
counts(:,:,i) = AllCompiledEmbryos{exp_idx(i)}.WindowedProfiles.HisRFP25CTime.Control.count(:,:,3);
ControlledProfiles(:,:,:,i) = AllCompiledEmbryos{exp_idx(i)}.WindowedProfiles.HisRFP25CTime.Control.mean;
end

FitInfo = cell(NumSets, NumSets, NChannels);
ScalingFactors = NaN(NumSets, NumSets, NChannels);
ScalingSEs = NaN(NumSets, NumSets, NChannels);
ScalingIntercepts = NaN(NumSets, NumSets, NChannels);
ScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
BestOverlapIndices = cell(NumSets, NumSets, NChannels);

GoodFitInfo = cell(NumSets, NumSets, NChannels);
GoodScalingFactors = NaN(NumSets, NumSets, NChannels);
GoodScalingSEs = NaN(NumSets, NumSets, NChannels);
GoodScalingIntercepts = NaN(NumSets, NumSets, NChannels);
GoodScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
GoodOverlapIndices = cell(NumSets, NumSets, NChannels);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
OkScalingIntercepts = NaN(NumSets, NumSets, NChannels);
OkScalingInterceptSEs = NaN(NumSets, NumSets, NChannels);
OkOverlapIndices = cell(NumSets, NumSets, NChannels);
for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            BestOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 5) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 5);
            GoodOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 3) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 3);
            OkOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 2) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 2);
            if sum(sum((OkOverlapIndices{exp2, exp1, ch_index} )) >= 1) >= 1
                
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1, ch_index}).';
        
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
               FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1, ch_index}).';
       
                dlm = fitlm(FitSet2, FitSet1);%,'Weights',  1./(FitSD1.^2));
                OkFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(2);
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(2);
                OkScalingIntercepts(exp2, exp1, ch_index) = dlm.Coefficients.Estimate(1);
                OkScalingInterceptSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE(1);

            end
        end
        
    end
end


for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            ScalingFactors(exp1, exp2, ch_index) = 1/ScalingFactors(exp2, exp1, ch_index);
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/(ScalingFactors(exp2, exp1, ch_index)^2);
            ScalingIntercepts(exp1, exp2, ch_index) = -ScalingIntercepts(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            ScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(ScalingInterceptSEs(exp2, exp1, ch_index)^2+ScalingIntercepts(exp1, exp2, ch_index)^2*ScalingSEs(exp2, exp1, ch_index)^2)/ScalingFactors(exp2, exp1, ch_index);
            GoodScalingFactors(exp1, exp2, ch_index) = 1/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingSEs(exp1, exp2, ch_index) = GoodScalingSEs(exp2, exp1, ch_index)/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingIntercepts(exp1, exp2, ch_index) = -GoodScalingIntercepts(exp2, exp1, ch_index)/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(GoodScalingInterceptSEs(exp2, exp1, ch_index)^2+GoodScalingIntercepts(exp1, exp2, ch_index)^2*GoodScalingSEs(exp2, exp1, ch_index)^2)/GoodScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingIntercepts(exp1, exp2, ch_index) = -OkScalingIntercepts(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingInterceptSEs(exp1, exp2, ch_index) = sqrt(OkScalingInterceptSEs(exp2, exp1, ch_index)^2+OkScalingIntercepts(exp1, exp2, ch_index)^2*OkScalingSEs(exp2, exp1, ch_index)^2)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingIntercepts = NaN(NumSets, NChannels);
SetScalingInterceptSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSet, 2:NChannels) = 1;
SetScalingSEs(RefSet, 2:NChannels) = 0;
SetScalingIntercepts(RefSet, 2:NChannels) = 0;
SetScalingInterceptSEs(RefSet, 2:NChannels) = 0;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, RefSet))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSet, ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSet, ch_index);
        SetScalingIntercepts(i, ch_index) = ScalingIntercepts(i, RefSet, ch_index);
        SetScalingInterceptSEs(i, ch_index) = ScalingInterceptSEs(i, RefSet, ch_index);
        if isnan(SetScalingFactors(i, ch_index) )
            SetScalingFactors(i, ch_index) = GoodScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = GoodScalingSEs(i, RefSet, ch_index);
            SetScalingIntercepts(i, ch_index) = GoodScalingIntercepts(i, RefSet, ch_index);
            SetScalingInterceptSEs(i, ch_index) = GoodScalingInterceptSEs(i, RefSet, ch_index);
        end
        if isnan(SetScalingFactors(i, ch_index) )
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSet, ch_index);
            SetScalingIntercepts(i, ch_index) = OkScalingIntercepts(i, RefSet, ch_index);
            SetScalingInterceptSEs(i, ch_index) = OkScalingInterceptSEs(i, RefSet, ch_index);
        end

    end
end

for exp_index = 1:length(AllSetInfo.Temperatures)
SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
CompiledEmbryos = AllCompiledEmbryos{exp_index};
CompiledEmbryos.WindowedHisRFP25CTimesControlScalingFactors = SetScalingFactors(exp_index,:);
CompiledEmbryos.WindowedHisRFP25CTimesControlScalingSEs = SetScalingSEs(exp_index,:);
CompiledEmbryos.WindowedHisRFP25CTimesControlScalingIntercepts = SetScalingIntercepts(exp_index,:);
CompiledEmbryos.WindowedHisRFP25CTimesControlScalingInterceptSEs = SetScalingInterceptSEs(exp_index,:);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
save(CEoutpath, 'CompiledEmbryos');
end

%%
% factor1 = [];
% factor2 = [];
% factor3 = [];
% channels = [];
% triple_multiplier = [];
% triple_multiplier_se = [];
% AllSets = 1:NumSets;
% for i = AllSets
%     for j = AllSets(~ismember(AllSets, [i]))
%         for k = AllSets(~ismember(AllSets, [i, j]))
%             for ch_index = 2:NChannels
%             factor1(end+1) = i;
%             factor2(end+1) = j;
%             factor3(end+1) = k;
%             channels(end+1) = ch_index;
%             triple_multiplier(end+1) = ScalingFactors(j, i, ch_index)*ScalingFactors(k, j, ch_index)*ScalingFactors(i, k, ch_index);
%             triple_multiplier_se(end+1) = sqrt(ScalingFactors(j, i, ch_index)^2*ScalingFactors(k, j, ch_index)^2*ScalingSEs(i, k, ch_index)^2 +...
%                 ScalingFactors(j, i, ch_index)^2*ScalingSEs(k, j, ch_index)^2*ScalingFactors(i, k, ch_index)^2+...
%                 ScalingSEs(j, i, ch_index)^2*ScalingFactors(k, j, ch_index)^2*ScalingFactors(i, k, ch_index)^2);
%             end
%         end
%     end
% end
% 
% % figure(2)
% % for i = AllSets
% %     scatter(i*ones(1, length(triple_multiplier(channels == 3 & (factor1 == i)))), triple_multiplier(channels == 3 & (factor1 == i)))
% %     hold on 
% % end
% % hold off
% % 
% % figure(3)
% % for i = AllSets
% %     scatter(i*ones(1, length(triple_multiplier(channels == 5 & (factor1 == i)))), triple_multiplier(channels == 5 & (factor1 == i)))
% %     hold on 
% % end
% % hold off
% 
% avg_triple_multiplier = NaN(NumSets,NChannels);
% bestsets = zeros(1, NChannels);
% for ch_index = 2:NChannels
%     for i = AllSets
%         
%         avg_triple_multiplier(i, ch_index) = mean(triple_multiplier(channels == ch_index & factor1 == i), 'omitnan');
%     end
%     [min_val, min_idx] = min(avg_triple_multiplier(:, ch_index));
%     bestsets(ch_index) = min_idx;
% end
% 
% SetScalingFactors = NaN(NumSets, NChannels);
% GoodFactors = zeros(NumSets, NumSets, NChannels, 'logical');
% % for ch_index = [3]%2:NChannels
% ch_index = 3;
% 
% BrokenCondition = false;
% IncludedSets = AllSets(any(~isnan(ScalingFactors(:,:,ch_index)), 1));
% SetsLeft = IncludedSets;
% SetsFound = [];
% 
% tm_notnan = ~isnan(triple_multiplier);
% factor1_sub = factor1(channels == ch_index & tm_notnan);
% factor2_sub = factor2(channels == ch_index & tm_notnan);
% factor3_sub = factor3(channels == ch_index & tm_notnan);
% 
% tm_sub = triple_multiplier(channels == ch_index & tm_notnan) ;
% tm_se_sub = triple_multiplier_se(channels == ch_index & tm_notnan) ;
% abs_tm_sub_diff = abs(tm_sub-1);
% [~, min_idx] = min(abs_tm_sub_diff);
% f1 = factor1_sub(min_idx);
% f2 = factor2_sub(min_idx);
% f3 = factor3_sub(min_idx);
% GoodFactors([f1 f1 f2 f2 f3 f3], [f2 f3 f3 f1 f1 f2], ch_index) = true;
% SetsFound = [SetsFound f1 f2 f3];
% SetsLeft = SetsLeft(~ismember(SetsLeft, SetsFound));
% RefSet = f1;
% SetScalingFactors(RefSet, ch_index) = 1;
% 
% counter = 0;
% while ~isempty(SetsLeft) & BrokenCondition == false
%     rm_cond = ((factor1_sub == f1)  & (factor2_sub == f2) & (factor3_sub == f3) ) | ...
%     ((factor1_sub == f1)  & (factor2_sub == f3) & (factor3_sub == f2) ) | ...
%     ((factor1_sub == f2)  & (factor2_sub == f3) & (factor3_sub == f1) ) | ...
%     ((factor1_sub == f2)  & (factor2_sub == f1) & (factor3_sub == f3) ) | ...
%     ((factor1_sub == f3)  & (factor2_sub == f1) & (factor3_sub == f2) ) | ...
%     ((factor1_sub == f3)  & (factor2_sub == f2) & (factor3_sub == f1) ) ;
% factor1_sub = factor1_sub(~rm_cond);
% factor2_sub = factor2_sub(~rm_cond);
% factor3_sub = factor3_sub(~rm_cond);
% 
% tm_sub = tm_sub(~rm_cond);
% tm_se_sub = tm_se_sub(~rm_cond);
% 
% is_good_tm = ismember(factor1_sub, SetsFound) | ismember(factor2_sub, SetsFound) | ismember(factor3_sub, SetsFound);
% if sum(is_good_tm) > 0
% factor1_temp = factor1_sub(is_good_tm);
% factor2_temp = factor2_sub(is_good_tm);
% factor3_temp = factor3_sub(is_good_tm);
% 
% tm_temp = tm_sub(is_good_tm);
% tm_se_temp = tm_se_sub(is_good_tm);
% 
% 
% abs_tm_temp_diff = abs(tm_temp-1);
% [min_val, min_idx] = min(abs_tm_temp_diff);
% f1 = factor1_temp(min_idx);
% f2 = factor2_temp(min_idx);
% f3 = factor3_temp(min_idx);
% GoodFactors(f1, f2, ch_index) = true;
% GoodFactors(f1, f3, ch_index) = true;
% GoodFactors(f2, f3, ch_index) = true;
% GoodFactors(f2, f1, ch_index) = true;
% GoodFactors(f3, f1, ch_index) = true;
% GoodFactors(f3, f2, ch_index) = true;
% NewFoundSets = [f1 f2 f3];
% SetsFound = [SetsFound NewFoundSets(~ismember(NewFoundSets, SetsFound))];
% SetsLeft = SetsLeft(~ismember(SetsLeft, SetsFound));
% else
%     BrokenCondition = true;
% end
% counter = counter + 1;
% end
% 
% 
