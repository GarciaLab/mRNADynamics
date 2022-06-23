function AllCompiledEmbryos = CalculateSetRescalingFactors
%%
SizeDataPath = 'S:/Gabriella/Dropbox/EmbryoSizeMeasurements/EmbryoSizeData.mat';
AllSetsProfFigPath = 'S:/Gabriella/Dropbox/ProteinProfiles/Figures/';
AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';
load(SizeDataPath,'APLengths','DVLengths','VentralDistances','DorsalDistances','NGoodEmbryos',...
    'MeanAPLength','StdAPLength','MeanDVLength','StdDVLength','AspectRatios');
MeanAspectRatio = mean(AspectRatios);
StdAspectRatio = std(AspectRatios);
[DubuisTimes, DubuisMeanProfile, ~] = getMembraneFurrowProfiles( 'dubuis');
[yw25CTimes, yw25CProfile, yw25CSE] = getMembraneFurrowProfiles( 'yw25csquished_nopv');
[hisrfp25CTimes, hisrfp25CProfile, hisrfp25CSE] = getMembraneFurrowProfiles( 'hisrfp25c_nopv');
AllSetInfo = GetFixedSetPrefixInfo;
 
NumSets = length(AllSetInfo.Temperatures);

NChannels = 5;


RefSets = [14, 14, 14, 14];
APbins = 0:0.025:1;
AllCompiledEmbryos = {};

ControlFittingBins = zeros(NChannels, length(APbins), 'logical');
ControlFittingBins(2,:) = APbins >= 0.2 & APbins <= 0.8; % Hoechst
ControlFittingBins(3,:) = APbins >= 0.1 & APbins <= 0.9 ;%APbins >= 0.1 & APbins <= 0.5 ; % Bicoid
ControlFittingBins(4,:) = APbins >= 0.1 & APbins <= 0.9 ;  %APbins >= 0.5 & APbins <= 0.9 ;  % Knirps
ControlFittingBins(5,:) = APbins >= 0.1 & APbins <= 0.9 ;% (APbins >= 0.1 & APbins <= 0.5) | (APbins >= 0.7 & APbins <= 0.9); % Hunchback
%%
for exp_index = 1:length(AllSetInfo.Temperatures)
if AllSetInfo.Flipped(exp_index)
    FlipString = 'yes';
else
    FlipString = 'no';
end
disp(['T = ', num2str(AllSetInfo.Temperatures(exp_index)), ', Rep: ', num2str(AllSetInfo.Replicates(exp_index)), ', Flipped: ', FlipString])
SetLabel = AllSetInfo.SetLabels{exp_index};
PlotLabel = AllSetInfo.PlotLabels{exp_index};
SetPrefixes = AllSetInfo.Prefixes{exp_index};
SetIsFlipped = AllSetInfo.Flipped(exp_index);
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
FixedPixelSize_um = liveExperiments{1}.pixelSize_um;
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
load(CEoutpath, 'CompiledEmbryos');
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
end




%%
warning('off','stats:LinearModel:RankDefDesignMat')
exp_idx = 1:NumSets;
AllSets = 1:NumSets;
counts_1sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_1sigma, 2));
counts_2sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_2sigma, 2));
SmoothedDeltaFCs = AllCompiledEmbryos{1}.FixCorrectedSmoothedDeltaFCs;
ControlledProfiles = NaN([size(AllCompiledEmbryos{1}.FixCorrectedSmoothedAvgAPProfiles.Control), NumSets]);
for i = 1:NumSets
counts_1sigma(i,:) = AllCompiledEmbryos{exp_idx(i)}.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_1sigma;
counts_2sigma(i,:) = AllCompiledEmbryos{exp_idx(i)}.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_2sigma;
ControlledProfiles(:,:,:,i) = AllCompiledEmbryos{exp_idx(i)}.FixCorrectedSmoothedAvgAPProfiles.Control;
end
FitInfo = cell(NumSets, NumSets, NChannels);
ScalingFactors = NaN(NumSets, NumSets, NChannels);
ScalingSEs = NaN(NumSets, NumSets, NChannels);
GoodOverlapIndices = cell(NumSets, NumSets);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
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
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                FitInfo{exp2, exp1, ch_index} = dlm;
                ScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                ScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
            end
        else
            disp([num2str(exp1),', ', num2str(exp2)])
        end
        
        if sum(OkOverlapIndices{exp2, exp1}) >= 8
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1},:);
                FitSet1 = FitSet1(:);
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
                FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1},:);
                FitSet2 = FitSet2(:);
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                OkFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
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
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSets(1), 2:5) = 1;
SetScalingSEs(RefSets(1), 2:5) = 1;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, [RefSets(1)]))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSets(1), ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSets(1), ch_index);
        if isnan(ScalingFactors(i, RefSets(1), ch_index))
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSets(1), ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSets(1), ch_index);
        end
    end
end

for exp_index = 1:length(AllSetInfo.Temperatures)
disp(['T = ', num2str(AllSetInfo.Temperatures(exp_index)), ', Rep: ', num2str(AllSetInfo.Replicates(exp_index)), ', Flipped: ', FlipString])
SetLabel = AllSetInfo.SetLabels{exp_index};
PlotLabel = AllSetInfo.PlotLabels{exp_index};
SetPrefixes = AllSetInfo.Prefixes{exp_index};
SetIsFlipped = AllSetInfo.Flipped(exp_index);
ProfFigPath = [AllSetsProfFigPath, SetLabel];
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
load(CEoutpath, 'CompiledEmbryos');
CompiledEmbryos.FixCorrectedControlScalingFactors = SetScalingFactors(exp_index,:);
CompiledEmbryos.FixCorrectedControlScalingSEs = SetScalingSEs(exp_index,:);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
save(CEoutpath, 'CompiledEmbryos');
end

%% Fit Universal Scaling using Dubuis Embryo Times
exp_idx = 1:NumSets;
AllSets = 1:NumSets;
counts_1sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_1sigma, 2));
counts_2sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_2sigma, 2));
counts_3sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_2sigma, 2));
SmoothedDeltaFCs = AllCompiledEmbryos{1}.DubuisSmoothedTimes;
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
GoodOverlapIndices = cell(NumSets, NumSets);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
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
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                FitInfo{exp2, exp1, ch_index} = dlm;
                ScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                ScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
            end
        else
            disp([num2str(exp1),', ', num2str(exp2)])
        end
        
        if sum(OkOverlapIndices{exp2, exp1}) >= 8
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1},:);
                FitSet1 = FitSet1(:);
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
                FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1},:);
                FitSet2 = FitSet2(:);
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                OkFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
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
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSets(2), 2:5) = 1;
SetScalingSEs(RefSets(2), 2:5) = 1;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, [RefSets(2)]))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSets(2), ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSets(2), ch_index);
        if isnan(ScalingFactors(i, RefSets(2), ch_index))
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSets(2), ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSets(2), ch_index);
        end
    end
end

for exp_index = 1:length(AllSetInfo.Temperatures)
if AllSetInfo.Flipped(exp_index)
    FlipString = 'yes';
else
    FlipString = 'no';
end
disp(['T = ', num2str(AllSetInfo.Temperatures(exp_index)), ', Rep: ', num2str(AllSetInfo.Replicates(exp_index)), ', Flipped: ', FlipString])
SetLabel = AllSetInfo.SetLabels{exp_index};
PlotLabel = AllSetInfo.PlotLabels{exp_index};
SetPrefixes = AllSetInfo.Prefixes{exp_index};
SetIsFlipped = AllSetInfo.Flipped(exp_index);
ProfFigPath = [AllSetsProfFigPath, SetLabel];
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
load(CEoutpath, 'CompiledEmbryos');
CompiledEmbryos.DubuisTimesControlScalingFactors = SetScalingFactors(exp_index,:);
CompiledEmbryos.DubuisTimesControlScalingSEs = SetScalingSEs(exp_index,:);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
save(CEoutpath, 'CompiledEmbryos');
end


%% Fit Universal Scaling using HisRFP 25C Times
exp_idx = 1:NumSets;
AllSets = 1:NumSets;
counts_1sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.HisRFPTimesSmoothedProfiles.ControlCounts.counts_within_1sigma, 2));
counts_2sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.HisRFPTimesSmoothedProfiles.ControlCounts.counts_within_2sigma, 2));
counts_3sigma = zeros(NumSets, size(AllCompiledEmbryos{1}.HisRFPTimesSmoothedProfiles.ControlCounts.counts_within_2sigma, 2));
SmoothedDeltaFCs = AllCompiledEmbryos{1}.HisRFP25CSmoothedTimes;
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
GoodOverlapIndices = cell(NumSets, NumSets);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
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
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                FitInfo{exp2, exp1, ch_index} = dlm;
                ScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                ScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
            end
        else
            disp([num2str(exp1),', ', num2str(exp2)])
        end
        
        if sum(OkOverlapIndices{exp2, exp1}) >= 8
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1},:);
                FitSet1 = FitSet1(:);
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
                FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1},:);
                FitSet2 = FitSet2(:);
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                OkFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
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
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSets(3), 2:5) = 1;
SetScalingSEs(RefSets(3), 2:5) = 1;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, [RefSets(3)]))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSets(3), ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSets(3), ch_index);
        if isnan(ScalingFactors(i, RefSets(3), ch_index))
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSets(3), ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSets(3), ch_index);
        end
    end
end

for exp_index = 1:length(AllSetInfo.Temperatures)
if AllSetInfo.Flipped(exp_index)
    FlipString = 'yes';
else
    FlipString = 'no';
end
disp(['T = ', num2str(AllSetInfo.Temperatures(exp_index)), ', Rep: ', num2str(AllSetInfo.Replicates(exp_index)), ', Flipped: ', FlipString])
SetLabel = AllSetInfo.SetLabels{exp_index};
PlotLabel = AllSetInfo.PlotLabels{exp_index};
SetPrefixes = AllSetInfo.Prefixes{exp_index};
SetIsFlipped = AllSetInfo.Flipped(exp_index);
ProfFigPath = [AllSetsProfFigPath, SetLabel];
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
load(CEoutpath, 'CompiledEmbryos');
CompiledEmbryos.HisRFP25CTimesControlScalingFactors = SetScalingFactors(exp_index,:);
CompiledEmbryos.HisRFP25CTimesControlScalingSEs = SetScalingSEs(exp_index,:);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
save(CEoutpath, 'CompiledEmbryos');
end

%% Windowed Dubuis Time Profile Info 

warning('off','stats:LinearModel:RankDefDesignMat')
RefSet = 14;
exp_idx = 1:NumSets;
AllSets = 1:NumSets;
WindowedDubuisTimes = AllCompiledEmbryos{1}.WindowedProfiles.DubuisTime.x;
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
BestOverlapIndices = cell(NumSets, NumSets, NChannels);

GoodFitInfo = cell(NumSets, NumSets, NChannels);
GoodScalingFactors = NaN(NumSets, NumSets, NChannels);
GoodScalingSEs = NaN(NumSets, NumSets, NChannels);
GoodOverlapIndices = cell(NumSets, NumSets, NChannels);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
OkOverlapIndices = cell(NumSets, NumSets, NChannels);
for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            BestOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 5) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 5);
            GoodOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 3) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 3);
            OkOverlapIndices{exp2, exp1, ch_index} = (counts(:,ControlFittingBins(ch_index, :),i) >= 2) &  (counts(:,ControlFittingBins(ch_index, :),j) >= 2);
            if sum(sum((BestOverlapIndices{exp2, exp1, ch_index} )) >= 5) >= 1
                
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( BestOverlapIndices{exp2, exp1, ch_index}).';
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
                FitSet2 = FitSet2( BestOverlapIndices{exp2, exp1, ch_index}).';
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                FitInfo{exp2, exp1, ch_index} = dlm;
                ScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                ScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
                
            end
            if sum(sum((GoodOverlapIndices{exp2, exp1, ch_index} )) >= 5) >= 1
                
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( GoodOverlapIndices{exp2, exp1, ch_index}).';
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
                FitSet2 = FitSet2( GoodOverlapIndices{exp2, exp1, ch_index}).';
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                GoodFitInfo{exp2, exp1, ch_index} = dlm;
                GoodScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                GoodScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
                
            end
            if sum(sum((OkOverlapIndices{exp2, exp1, ch_index} )) >= 2) >= 1
                
                FitSet1 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1, ch_index}).';
                FitSet2 = ControlledProfiles(:,ControlFittingBins(ch_index, :),ch_index,j);
                FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1, ch_index}).';
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                OkFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
                
            else
                disp([num2str(exp1),', ', num2str(exp2)])
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
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            GoodScalingFactors(exp1, exp2, ch_index) = 1/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingSEs(exp1, exp2, ch_index) = GoodScalingSEs(exp2, exp1, ch_index)/GoodScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSet, 2:5) = 1;
SetScalingSEs(RefSet, 2:5) = 1;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, [RefSet]))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSet, ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSet, ch_index);
        if isnan(ScalingFactors(i, RefSet, ch_index))
            SetScalingFactors(i, ch_index) = GoodScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = GoodScalingSEs(i, RefSet, ch_index);
        end
        
        if isnan(ScalingFactors(i, RefSet, ch_index))
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSet, ch_index);
        end
    end
end

for exp_index = 1:length(AllSetInfo.Temperatures)
disp(['T = ', num2str(AllSetInfo.Temperatures(exp_index)), ', Rep: ', num2str(AllSetInfo.Replicates(exp_index)), ', Flipped: ', FlipString])
SetLabel = AllSetInfo.SetLabels{exp_index};
PlotLabel = AllSetInfo.PlotLabels{exp_index};
SetPrefixes = AllSetInfo.Prefixes{exp_index};
SetIsFlipped = AllSetInfo.Flipped(exp_index);
ProfFigPath = [AllSetsProfFigPath, SetLabel];
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
load(CEoutpath, 'CompiledEmbryos');
CompiledEmbryos.WindowedDubuisTimesControlScalingFactors = SetScalingFactors(exp_index,:);
CompiledEmbryos.WindowedDubuisTimesControlScalingSEs = SetScalingSEs(exp_index,:);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
save(CEoutpath, 'CompiledEmbryos');
end

%% Windowed Delta FC Profile Info 

warning('off','stats:LinearModel:RankDefDesignMat')
RefSet = 14;
exp_idx = 1:NumSets;
AllSets = 1:NumSets;
WindowedDeltaFCTimes = AllCompiledEmbryos{1}.WindowedProfiles.DeltaFC.x;
CountDims = [size(AllCompiledEmbryos{1}.WindowedProfiles.DeltaFC.Control.count, 1),...
    size(AllCompiledEmbryos{1}.WindowedProfiles.DeltaFC.Control.count, 2), NumSets];
counts = zeros(CountDims);

ControlledProfiles = NaN([size(AllCompiledEmbryos{1}.WindowedProfiles.DeltaFC.Control.count) NumSets]);
for i = 1:NumSets
counts(:,:,i) = AllCompiledEmbryos{exp_idx(i)}.WindowedProfiles.DeltaFC.Control.count(:,:,3);
ControlledProfiles(:,:,:,i) = AllCompiledEmbryos{exp_idx(i)}.WindowedProfiles.DeltaFC.Control.mean;
end
FitInfo = cell(NumSets, NumSets, NChannels);
ScalingFactors = NaN(NumSets, NumSets, NChannels);
ScalingSEs = NaN(NumSets, NumSets, NChannels);
BestOverlapIndices = cell(NumSets, NumSets);

GoodFitInfo = cell(NumSets, NumSets, NChannels);
GoodScalingFactors = NaN(NumSets, NumSets, NChannels);
GoodScalingSEs = NaN(NumSets, NumSets, NChannels);
GoodOverlapIndices = cell(NumSets, NumSets);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
OkOverlapIndices = cell(NumSets, NumSets);
for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        BestOverlapIndices{exp2, exp1} = (counts(:,:,i) >= 5) &  (counts(:,:,j) >= 5);
        GoodOverlapIndices{exp2, exp1} = (counts(:,:,i) >= 3) &  (counts(:,:,j) >= 3);
        OkOverlapIndices{exp2, exp1} = (counts(:,:,i) >= 2) &  (counts(:,:,j) >= 2);
        if sum(sum((BestOverlapIndices{exp2, exp1} )) >= 5) >= 20
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,:,ch_index,i);
                FitSet1 = FitSet1( BestOverlapIndices{exp2, exp1}).';
                FitSet2 = ControlledProfiles(:,:,ch_index,j);
                FitSet2 = FitSet2( BestOverlapIndices{exp2, exp1}).';
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                FitInfo{exp2, exp1, ch_index} = dlm;
                ScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                ScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
            end
        end
        if sum(sum((GoodOverlapIndices{exp2, exp1} )) >= 5) >= 20
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,:,ch_index,i);
                FitSet1 = FitSet1( GoodOverlapIndices{exp2, exp1}).';
                FitSet2 = ControlledProfiles(:,:,ch_index,j);
                FitSet2 = FitSet2( GoodOverlapIndices{exp2, exp1}).';
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                GoodFitInfo{exp2, exp1, ch_index} = dlm;
                GoodScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                GoodScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
            end
        end
        if sum(sum((OkOverlapIndices{exp2, exp1} )) >= 2) >= 20
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,:,ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1}).';
                FitSet2 = ControlledProfiles(:,:,ch_index,j);
                FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1}).';
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                OkGoodFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
            end
        else
            disp([num2str(exp1),', ', num2str(exp2)])
        end
        
     
    end
end

for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            ScalingFactors(exp1, exp2, ch_index) = 1/ScalingFactors(exp2, exp1, ch_index);
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            GoodScalingFactors(exp1, exp2, ch_index) = 1/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingSEs(exp1, exp2, ch_index) = GoodScalingSEs(exp2, exp1, ch_index)/GoodScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSet, 2:5) = 1;
SetScalingSEs(RefSet, 2:5) = 1;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, [RefSet]))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSet, ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSet, ch_index);
        if isnan(ScalingFactors(i, RefSet, ch_index))
            SetScalingFactors(i, ch_index) = GoodScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = GoodScalingSEs(i, RefSet, ch_index);
        end
        
        if isnan(ScalingFactors(i, RefSet, ch_index))
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSet, ch_index);
        end
    end
end

for exp_index = 1:length(AllSetInfo.Temperatures)
disp(['T = ', num2str(AllSetInfo.Temperatures(exp_index)), ', Rep: ', num2str(AllSetInfo.Replicates(exp_index)), ', Flipped: ', FlipString])
SetLabel = AllSetInfo.SetLabels{exp_index};
PlotLabel = AllSetInfo.PlotLabels{exp_index};
SetPrefixes = AllSetInfo.Prefixes{exp_index};
SetIsFlipped = AllSetInfo.Flipped(exp_index);
ProfFigPath = [AllSetsProfFigPath, SetLabel];
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
load(CEoutpath, 'CompiledEmbryos');
CompiledEmbryos.WindowedDeltaFCControlScalingFactors = SetScalingFactors(exp_index,:);
CompiledEmbryos.WindowedDeltaFCControlScalingSEs = SetScalingSEs(exp_index,:);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
save(CEoutpath, 'CompiledEmbryos');
end

%% Windowed yw25C Profile Info 

warning('off','stats:LinearModel:RankDefDesignMat')
RefSet = 14;
exp_idx = 1:NumSets;
AllSets = 1:NumSets;
WindowedDeltaFCTimes = AllCompiledEmbryos{1}.WindowedProfiles.yw25CTime.x;
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
BestOverlapIndices = cell(NumSets, NumSets);

GoodFitInfo = cell(NumSets, NumSets, NChannels);
GoodScalingFactors = NaN(NumSets, NumSets, NChannels);
GoodScalingSEs = NaN(NumSets, NumSets, NChannels);
GoodOverlapIndices = cell(NumSets, NumSets);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
OkOverlapIndices = cell(NumSets, NumSets);
for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        BestOverlapIndices{exp2, exp1} = (counts(:,:,i) >= 5) &  (counts(:,:,j) >= 5);
        GoodOverlapIndices{exp2, exp1} = (counts(:,:,i) >= 3) &  (counts(:,:,j) >= 3);
        OkOverlapIndices{exp2, exp1} = (counts(:,:,i) >= 2) &  (counts(:,:,j) >= 2);
        if sum(sum((BestOverlapIndices{exp2, exp1} )) >= 5) >= 20
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,:,ch_index,i);
                FitSet1 = FitSet1( BestOverlapIndices{exp2, exp1}).';
                FitSet2 = ControlledProfiles(:,:,ch_index,j);
                FitSet2 = FitSet2( BestOverlapIndices{exp2, exp1}).';
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                FitInfo{exp2, exp1, ch_index} = dlm;
                ScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                ScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
            end
        end
        if sum(sum((GoodOverlapIndices{exp2, exp1} )) >= 5) >= 20
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,:,ch_index,i);
                FitSet1 = FitSet1( GoodOverlapIndices{exp2, exp1}).';
                FitSet2 = ControlledProfiles(:,:,ch_index,j);
                FitSet2 = FitSet2( GoodOverlapIndices{exp2, exp1}).';
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                GoodFitInfo{exp2, exp1, ch_index} = dlm;
                GoodScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                GoodScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
            end
        end
        if sum(sum((OkOverlapIndices{exp2, exp1} )) >= 2) >= 20
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,:,ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1}).';
                FitSet2 = ControlledProfiles(:,:,ch_index,j);
                FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1}).';
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                OkGoodFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
            end
        else
            disp([num2str(exp1),', ', num2str(exp2)])
        end
        
     
    end
end

for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            ScalingFactors(exp1, exp2, ch_index) = 1/ScalingFactors(exp2, exp1, ch_index);
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            GoodScalingFactors(exp1, exp2, ch_index) = 1/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingSEs(exp1, exp2, ch_index) = GoodScalingSEs(exp2, exp1, ch_index)/GoodScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSet, 2:5) = 1;
SetScalingSEs(RefSet, 2:5) = 1;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, [RefSet]))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSet, ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSet, ch_index);
        if isnan(ScalingFactors(i, RefSet, ch_index))
            SetScalingFactors(i, ch_index) = GoodScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = GoodScalingSEs(i, RefSet, ch_index);
        end
        
        if isnan(ScalingFactors(i, RefSet, ch_index))
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSet, ch_index);
        end
    end
end

for exp_index = 1:length(AllSetInfo.Temperatures)
disp(['T = ', num2str(AllSetInfo.Temperatures(exp_index)), ', Rep: ', num2str(AllSetInfo.Replicates(exp_index)), ', Flipped: ', FlipString])
SetLabel = AllSetInfo.SetLabels{exp_index};
PlotLabel = AllSetInfo.PlotLabels{exp_index};
SetPrefixes = AllSetInfo.Prefixes{exp_index};
SetIsFlipped = AllSetInfo.Flipped(exp_index);
ProfFigPath = [AllSetsProfFigPath, SetLabel];
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
load(CEoutpath, 'CompiledEmbryos');
CompiledEmbryos.Windowedyw25CTimesControlScalingFactors = SetScalingFactors(exp_index,:);
CompiledEmbryos.Windowedyw25CTimesControlScalingSEs = SetScalingSEs(exp_index,:);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
save(CEoutpath, 'CompiledEmbryos');
end

%% Windowed HisRFP 25C Profile Info 

warning('off','stats:LinearModel:RankDefDesignMat')
RefSet = 14;
exp_idx = 1:NumSets;
AllSets = 1:NumSets;
WindowedDeltaFCTimes = AllCompiledEmbryos{1}.WindowedProfiles.HisRFP25CTime.x;
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
BestOverlapIndices = cell(NumSets, NumSets);

GoodFitInfo = cell(NumSets, NumSets, NChannels);
GoodScalingFactors = NaN(NumSets, NumSets, NChannels);
GoodScalingSEs = NaN(NumSets, NumSets, NChannels);
GoodOverlapIndices = cell(NumSets, NumSets);

OkFitInfo = cell(NumSets, NumSets, NChannels);
OkScalingFactors = NaN(NumSets, NumSets, NChannels);
OkScalingSEs = NaN(NumSets, NumSets, NChannels);
OkOverlapIndices = cell(NumSets, NumSets);
for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        BestOverlapIndices{exp2, exp1} = (counts(:,:,i) >= 5) &  (counts(:,:,j) >= 5);
        GoodOverlapIndices{exp2, exp1} = (counts(:,:,i) >= 3) &  (counts(:,:,j) >= 3);
        OkOverlapIndices{exp2, exp1} = (counts(:,:,i) >= 2) &  (counts(:,:,j) >= 2);
        if sum(sum((BestOverlapIndices{exp2, exp1} )) >= 5) >= 20
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,:,ch_index,i);
                FitSet1 = FitSet1( BestOverlapIndices{exp2, exp1}).';
                FitSet2 = ControlledProfiles(:,:,ch_index,j);
                FitSet2 = FitSet2( BestOverlapIndices{exp2, exp1}).';
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                FitInfo{exp2, exp1, ch_index} = dlm;
                ScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                ScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
            end
        end
        if sum(sum((GoodOverlapIndices{exp2, exp1} )) >= 5) >= 20
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,:,ch_index,i);
                FitSet1 = FitSet1( GoodOverlapIndices{exp2, exp1}).';
                FitSet2 = ControlledProfiles(:,:,ch_index,j);
                FitSet2 = FitSet2( GoodOverlapIndices{exp2, exp1}).';
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                GoodFitInfo{exp2, exp1, ch_index} = dlm;
                GoodScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                GoodScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
            end
        end
        if sum(sum((OkOverlapIndices{exp2, exp1} )) >= 2) >= 20
            for ch_index = 2:NChannels
                FitSet1 = ControlledProfiles(:,:,ch_index,i);
                FitSet1 = FitSet1( OkOverlapIndices{exp2, exp1}).';
                FitSet2 = ControlledProfiles(:,:,ch_index,j);
                FitSet2 = FitSet2( OkOverlapIndices{exp2, exp1}).';
                dlm = fitlm(FitSet2, FitSet1, 'Intercept', false);
                OkGoodFitInfo{exp2, exp1, ch_index} = dlm;
                OkScalingFactors(exp2, exp1, ch_index) = dlm.Coefficients.Estimate;
                OkScalingSEs(exp2, exp1, ch_index) = dlm.Coefficients.SE;
            end
        else
            disp([num2str(exp1),', ', num2str(exp2)])
        end
        
     
    end
end

for i = 1:(NumSets-1)
    exp1 = exp_idx(i);
    for j = i+1:NumSets
        exp2 = exp_idx(j);
        for ch_index = 2:NChannels
            ScalingFactors(exp1, exp2, ch_index) = 1/ScalingFactors(exp2, exp1, ch_index);
            ScalingSEs(exp1, exp2, ch_index) = ScalingSEs(exp2, exp1, ch_index)/ScalingFactors(exp2, exp1, ch_index);
            GoodScalingFactors(exp1, exp2, ch_index) = 1/GoodScalingFactors(exp2, exp1, ch_index);
            GoodScalingSEs(exp1, exp2, ch_index) = GoodScalingSEs(exp2, exp1, ch_index)/GoodScalingFactors(exp2, exp1, ch_index);
            OkScalingFactors(exp1, exp2, ch_index) = 1/OkScalingFactors(exp2, exp1, ch_index);
            OkScalingSEs(exp1, exp2, ch_index) = OkScalingSEs(exp2, exp1, ch_index)/OkScalingFactors(exp2, exp1, ch_index);
        end
        
    end
end

SetScalingFactors = NaN(NumSets, NChannels);
SetScalingSEs = NaN(NumSets, NChannels);
SetScalingFactors(RefSet, 2:5) = 1;
SetScalingSEs(RefSet, 2:5) = 1;
for ch_index = 2:NChannels 
    for i = AllSets(~ismember(AllSets, [RefSet]))
        SetScalingFactors(i, ch_index) = ScalingFactors(i, RefSet, ch_index);
        SetScalingSEs(i, ch_index) = ScalingSEs(i, RefSet, ch_index);
        if isnan(ScalingFactors(i, RefSet, ch_index))
            SetScalingFactors(i, ch_index) = GoodScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = GoodScalingSEs(i, RefSet, ch_index);
        end
        
        if isnan(ScalingFactors(i, RefSet, ch_index))
            SetScalingFactors(i, ch_index) = OkScalingFactors(i, RefSet, ch_index);
            SetScalingSEs(i, ch_index) = OkScalingSEs(i, RefSet, ch_index);
        end
    end
end

for exp_index = 1:length(AllSetInfo.Temperatures)
disp(['T = ', num2str(AllSetInfo.Temperatures(exp_index)), ', Rep: ', num2str(AllSetInfo.Replicates(exp_index)), ', Flipped: ', FlipString])
SetLabel = AllSetInfo.SetLabels{exp_index};
PlotLabel = AllSetInfo.PlotLabels{exp_index};
SetPrefixes = AllSetInfo.Prefixes{exp_index};
SetIsFlipped = AllSetInfo.Flipped(exp_index);
ProfFigPath = [AllSetsProfFigPath, SetLabel];
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
load(CEoutpath, 'CompiledEmbryos');
CompiledEmbryos.WindowedHisRFP25CTimesControlScalingFactors = SetScalingFactors(exp_index,:);
CompiledEmbryos.WindowedHisRFP25CTimesControlScalingSEs = SetScalingSEs(exp_index,:);
AllCompiledEmbryos{exp_index} = CompiledEmbryos;
save(CEoutpath, 'CompiledEmbryos');
end
%%
factor1 = [];
factor2 = [];
factor3 = [];
channels = [];
triple_multiplier = [];
triple_multiplier_se = [];
AllSets = 1:NumSets;
for i = AllSets
    for j = AllSets(~ismember(AllSets, [i]))
        for k = AllSets(~ismember(AllSets, [i, j]))
            for ch_index = 2:NChannels
            factor1(end+1) = i;
            factor2(end+1) = j;
            factor3(end+1) = k;
            channels(end+1) = ch_index;
            triple_multiplier(end+1) = ScalingFactors(j, i, ch_index)*ScalingFactors(k, j, ch_index)*ScalingFactors(i, k, ch_index);
            triple_multiplier_se(end+1) = sqrt(ScalingFactors(j, i, ch_index)^2*ScalingFactors(k, j, ch_index)^2*ScalingSEs(i, k, ch_index)^2 +...
                ScalingFactors(j, i, ch_index)^2*ScalingSEs(k, j, ch_index)^2*ScalingFactors(i, k, ch_index)^2+...
                ScalingSEs(j, i, ch_index)^2*ScalingFactors(k, j, ch_index)^2*ScalingFactors(i, k, ch_index)^2);
            end
        end
    end
end

% figure(2)
% for i = AllSets
%     scatter(i*ones(1, length(triple_multiplier(channels == 3 & (factor1 == i)))), triple_multiplier(channels == 3 & (factor1 == i)))
%     hold on 
% end
% hold off
% 
% figure(3)
% for i = AllSets
%     scatter(i*ones(1, length(triple_multiplier(channels == 5 & (factor1 == i)))), triple_multiplier(channels == 5 & (factor1 == i)))
%     hold on 
% end
% hold off

avg_triple_multiplier = NaN(NumSets,NChannels);
bestsets = zeros(1, NChannels);
for ch_index = 2:NChannels
    for i = AllSets
        
        avg_triple_multiplier(i, ch_index) = mean(triple_multiplier(channels == ch_index & factor1 == i), 'omitnan');
    end
    [min_val, min_idx] = min(avg_triple_multiplier(:, ch_index));
    bestsets(ch_index) = min_idx;
end

SetScalingFactors = NaN(NumSets, NChannels);
GoodFactors = zeros(NumSets, NumSets, NChannels, 'logical');
% for ch_index = [3]%2:NChannels
ch_index = 3;

BrokenCondition = false;
IncludedSets = AllSets(any(~isnan(ScalingFactors(:,:,ch_index)), 1));
SetsLeft = IncludedSets;
SetsFound = [];

tm_notnan = ~isnan(triple_multiplier);
factor1_sub = factor1(channels == ch_index & tm_notnan);
factor2_sub = factor2(channels == ch_index & tm_notnan);
factor3_sub = factor3(channels == ch_index & tm_notnan);

tm_sub = triple_multiplier(channels == ch_index & tm_notnan) ;
tm_se_sub = triple_multiplier_se(channels == ch_index & tm_notnan) ;
abs_tm_sub_diff = abs(tm_sub-1);
[~, min_idx] = min(abs_tm_sub_diff);
f1 = factor1_sub(min_idx);
f2 = factor2_sub(min_idx);
f3 = factor3_sub(min_idx);
GoodFactors([f1 f1 f2 f2 f3 f3], [f2 f3 f3 f1 f1 f2], ch_index) = true;
SetsFound = [SetsFound f1 f2 f3];
SetsLeft = SetsLeft(~ismember(SetsLeft, SetsFound));
RefSet = f1;
SetScalingFactors(RefSet, ch_index) = 1;

counter = 0;
while ~isempty(SetsLeft) & BrokenCondition == false
    rm_cond = ((factor1_sub == f1)  & (factor2_sub == f2) & (factor3_sub == f3) ) | ...
    ((factor1_sub == f1)  & (factor2_sub == f3) & (factor3_sub == f2) ) | ...
    ((factor1_sub == f2)  & (factor2_sub == f3) & (factor3_sub == f1) ) | ...
    ((factor1_sub == f2)  & (factor2_sub == f1) & (factor3_sub == f3) ) | ...
    ((factor1_sub == f3)  & (factor2_sub == f1) & (factor3_sub == f2) ) | ...
    ((factor1_sub == f3)  & (factor2_sub == f2) & (factor3_sub == f1) ) ;
factor1_sub = factor1_sub(~rm_cond);
factor2_sub = factor2_sub(~rm_cond);
factor3_sub = factor3_sub(~rm_cond);

tm_sub = tm_sub(~rm_cond);
tm_se_sub = tm_se_sub(~rm_cond);

is_good_tm = ismember(factor1_sub, SetsFound) | ismember(factor2_sub, SetsFound) | ismember(factor3_sub, SetsFound);
if sum(is_good_tm) > 0
factor1_temp = factor1_sub(is_good_tm);
factor2_temp = factor2_sub(is_good_tm);
factor3_temp = factor3_sub(is_good_tm);

tm_temp = tm_sub(is_good_tm);
tm_se_temp = tm_se_sub(is_good_tm);


abs_tm_temp_diff = abs(tm_temp-1);
[min_val, min_idx] = min(abs_tm_temp_diff);
f1 = factor1_temp(min_idx);
f2 = factor2_temp(min_idx);
f3 = factor3_temp(min_idx);
GoodFactors(f1, f2, ch_index) = true;
GoodFactors(f1, f3, ch_index) = true;
GoodFactors(f2, f3, ch_index) = true;
GoodFactors(f2, f1, ch_index) = true;
GoodFactors(f3, f1, ch_index) = true;
GoodFactors(f3, f2, ch_index) = true;
NewFoundSets = [f1 f2 f3];
SetsFound = [SetsFound NewFoundSets(~ismember(NewFoundSets, SetsFound))];
SetsLeft = SetsLeft(~ismember(SetsLeft, SetsFound));
else
    BrokenCondition = true;
end
counter = counter + 1;
end


