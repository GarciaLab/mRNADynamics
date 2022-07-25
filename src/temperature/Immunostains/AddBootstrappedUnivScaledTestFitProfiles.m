function CompiledEmbryos =  AddBootstrappedUnivScaledTestFitProfiles(CompiledEmbryos, exp_index)
AllSetInfo = GetFixedSetPrefixInfo;

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

BootstrappedProfilePath = 'S:/Gabriella/Dropbox/BootstrappedTestData/25CBootstrappedGaussianSmoothedProfiles.mat';
load(BootstrappedProfilePath, 'MeanSmoothedProfiles', 'SmoothedProfileSEs', 'xfits');
NumMasterProfs = length(MeanSmoothedProfiles);

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];

[NumEmbryos, NumAPbins, NChannels] = size(CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles);
FitTypeStrings = {'Rep1Rep2ComboFitSlideRescaledDorsalAvgAPProfiles', 'Rep1FlippedComboFitSlideRescaledDorsalAvgAPProfiles',...
    'Rep2FlippedComboFitSlideRescaledDorsalAvgAPProfiles', 'AllCombinedFitSlideRescaledDorsalAvgAPProfiles',...
    'Rep1Rep2ComboFitZeroedSlideRescaledDorsalAvgAPProfiles', 'Rep1FlippedComboFitZeroedSlideRescaledDorsalAvgAPProfiles',...
    'Rep2FlippedComboFitZeroedSlideRescaledDorsalAvgAPProfiles', 'AllCombinedFitZeroedSlideRescaledDorsalAvgAPProfiles',...
    'Rep1Rep2ComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles', 'Rep1FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles',...
    'Rep2FlippedComboZeroCorrectedSlideRescaledDorsalAvgAPProfiles', 'AllCombinedZeroCorrectedSlideRescaledDorsalAvgAPProfiles',...
    'Rep1Rep2ComboSlideRescaledDorsalAvgAPProfiles', 'Rep1FlippedComboSlideRescaledDorsalAvgAPProfiles',...
    'Rep2FlippedComboSlideRescaledDorsalAvgAPProfiles', 'AllCombinedSlideRescaledDorsalAvgAPProfiles',...
    'Rep1Rep2ComboHybridRescaledDorsalAvgAPProfiles', 'Rep1FlippedComboHybridRescaledDorsalAvgAPProfiles',...
    'Rep2FlippedComboHybridRescaledDorsalAvgAPProfiles', 'AllCombinedHybridRescaledDorsalAvgAPProfiles'...
    };
SetStrings = {'Rep1', 'Rep2', 'Flipped'};
chLists = {[3], [3], [3], [3],...
    [3], [3], [3], [3],...
    [3 5], [3 5], [3 5], [3 5],...
    [3 5], [3 5], [3 5], [3 5],...
    [3], [3], [3], [3]};


%%

xfits = 0:1:70;
Nxfits = length(xfits);
min_2sigma_points = 5;
counts_above_limit = 3;
counts_below_limit = 3;
sigma = 5;
window_multiplier = 3;
NumBootstrappedFits = 200;
NumPoints = 100;


%% Bootstrapping the FitSlideRescaledDorsalAvgAPProfiles for the Test Set]
UseTestTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
UseControlTF = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);

if ~isfield(CompiledEmbryos, 'BootstrappedUnivScaledProfiles')
    CompiledEmbryos.BootstrappedUnivScaledProfiles = {};
end
for i = 1:length(FitTypeStrings)
    CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}) = {};
    for j = 1:length(SetStrings)
        disp(['i = ', num2str(i), '/', num2str(length(FitTypeStrings)),', j = ', num2str(j), '/', num2str(length(SetStrings))])
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}) = {};
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).x = xfits;
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).ControlSet = {};
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).ControlSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfs);
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).ControlSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfs);
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).ControlSet.counts = NaN(Nxfits, NumMasterProfs);
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).ControlSet.counts_above = NaN(Nxfits, NumMasterProfs);
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).ControlSet.counts_below = NaN(Nxfits, NumMasterProfs);
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).TestSet = {};
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).TestSet.mean = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfs);
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).TestSet.se = NaN(Nxfits, NumAPbins, NChannels, NumMasterProfs);
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).TestSet.counts = NaN(Nxfits, NumMasterProfs);
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).TestSet.counts_above = NaN(Nxfits, NumMasterProfs);
        CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).TestSet.counts_below = NaN(Nxfits, NumMasterProfs);
        for ch_index = chLists{i}
            for master_index = 1:NumMasterProfs
                x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTestTF);
                ys = CompiledEmbryos.UnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j});
                
                ys = ys(UseTestTF, :,ch_index, master_index);
                BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
                
                x_sample = x_sample(~BadTF);
                ys = ys(~BadTF,:);
                [x_sample, sort_order] = sort(x_sample);
                ys = ys(sort_order,:);
                
                
                DiffMat = xfits.'-x_sample;
                GaussianWeights = GetGaussianWeightMat(xfits, x_sample, sigma, window_multiplier*sigma);
                
                SmoothedProfiles = NaN(Nxfits, NumAPbins, NumBootstrappedFits);
                counts = zeros(1, Nxfits);
                counts_above = zeros(1,Nxfits);
                counts_below = zeros(1, Nxfits);
                
                for x_index = 1:Nxfits
                    MatchedPoints = find(GaussianWeights(x_index,:)> 0);
                    Num2SigmaPoints = sum(DiffMat(x_index,MatchedPoints) > -2*sigma & DiffMat(x_index,MatchedPoints) < 2*sigma );
                    
                    if Num2SigmaPoints < min_2sigma_points
                        GaussianWeights(x_index,:) = 0;
                        continue
                    end
                    
                    counts(x_index) = Num2SigmaPoints;
                    MatchedPoints = find(GaussianWeights(x_index,:)> 0);
                    
                    counts_above(x_index) = sum(DiffMat(x_index,MatchedPoints) > 0 & DiffMat(x_index,MatchedPoints) < 2*sigma );
                    counts_below(x_index) = sum(DiffMat(x_index,MatchedPoints) < 0 & DiffMat(x_index,MatchedPoints) > -2*sigma );
                    
                    if counts_above(x_index) < counts_above_limit
                        GaussianWeights(x_index,:) = 0;
                    end
                    
                    if counts_below(x_index) < counts_below_limit
                        GaussianWeights(x_index,:) = 0;
                    end
                    
                end
                
                
                ValidCountTFs = (counts >= min_2sigma_points & counts_above >= counts_above_limit & counts_below >= counts_below_limit & xfits >= min(x_sample)+sigma/2 & xfits <= max(x_sample) - sigma/2);
                x_indices = 1:Nxfits;
                ValidXfits = x_indices(ValidCountTFs);
                NValidXfits = length(ValidXfits);
                
                if NValidXfits > 0
                    for k = 1:NValidXfits
                        WindowWidthLimits = [xfits(ValidXfits(k))-( min(x_sample)), max(x_sample)- xfits(ValidXfits(k))];
                        for rep = 1:NumBootstrappedFits*2
                            r1  = (rand(NumPoints, 1)*(2*window_multiplier*sigma)-window_multiplier*sigma).';
                            Deltas = DiffMat(ValidXfits(k),:);
                            Deltas(Deltas > window_multiplier*sigma | Deltas < -window_multiplier*sigma) = NaN;
                            ind = NaN(1, NumPoints);
                            for l = 1:NumPoints
                                ind(l) = find(min(abs(Deltas - r1(l))) == abs(Deltas - r1(l)));
                            end
                            SmoothedProfiles(ValidXfits(k),:,rep) =  GaussianWeights(ValidXfits(k),ind)*ys(ind, :)/(sum(GaussianWeights(ValidXfits(k),ind)));
                        end
                    end
                    
                    
                    
                    MeanProf = mean(SmoothedProfiles, 3, 'omitnan');
                    SEProf = std(SmoothedProfiles,0,  3, 'omitnan')/sqrt(NumBootstrappedFits);
                    SEProf(MeanProf == 0) = NaN;
                    MeanProf(MeanProf == 0) = NaN;
                    
                    CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).x = xfits;
                    CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).TestSet.mean(:,:,ch_index, master_index) = MeanProf;
                    CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).TestSet.se(:,:,ch_index, master_index) = SEProf;
                    CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).TestSet.counts(:,master_index) = counts;
                    CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).TestSet.counts_above(:,master_index) = counts_above;
                    CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).TestSet.counts_below(:,master_index) = counts_below;
                end
                % Control Set Smoothing
                x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseControlTF);
                ys = CompiledEmbryos.UnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j});
                
                ys = ys(UseControlTF, :,ch_index, master_index);
                BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
                
                x_sample = x_sample(~BadTF);
                ys = ys(~BadTF,:);
                [x_sample, sort_order] = sort(x_sample);
                ys = ys(sort_order,:);
                
                
                DiffMat = xfits.'-x_sample;
                GaussianWeights = GetGaussianWeightMat(xfits, x_sample, sigma, window_multiplier*sigma);
                
                SmoothedProfiles = NaN(Nxfits, NumAPbins, NumBootstrappedFits);
                counts = zeros(1, Nxfits);
                counts_above = zeros(1,Nxfits);
                counts_below = zeros(1, Nxfits);
                
                for x_index = 1:Nxfits
                    MatchedPoints = find(GaussianWeights(x_index,:)> 0);
                    Num2SigmaPoints = sum(DiffMat(x_index,MatchedPoints) > -2*sigma & DiffMat(x_index,MatchedPoints) < 2*sigma );
                    
                    if Num2SigmaPoints < min_2sigma_points
                        GaussianWeights(x_index,:) = 0;
                        continue
                    end
                    
                    counts(x_index) = Num2SigmaPoints;
                    MatchedPoints = find(GaussianWeights(x_index,:)> 0);
                    
                    counts_above(x_index) = sum(DiffMat(x_index,MatchedPoints) > 0 & DiffMat(x_index,MatchedPoints) < 2*sigma );
                    counts_below(x_index) = sum(DiffMat(x_index,MatchedPoints) < 0 & DiffMat(x_index,MatchedPoints) > -2*sigma );
                    
                    if counts_above(x_index) < counts_above_limit
                        GaussianWeights(x_index,:) = 0;
                    end
                    
                    if counts_below(x_index) < counts_below_limit
                        GaussianWeights(x_index,:) = 0;
                    end
                    
                end
                
                
                ValidCountTFs = (counts >= min_2sigma_points & counts_above >= counts_above_limit & counts_below >= counts_below_limit & xfits >= min(x_sample)+sigma/2 & xfits <= max(x_sample) - sigma/2);
                x_indices = 1:Nxfits;
                ValidXfits = x_indices(ValidCountTFs);
                NValidXfits = length(ValidXfits);
                
                if NValidXfits > 0
                    for k = 1:NValidXfits
                        WindowWidthLimits = [xfits(ValidXfits(k))-( min(x_sample)), max(x_sample)- xfits(ValidXfits(k))];
                        for rep = 1:NumBootstrappedFits*2
                            r1  = (rand(NumPoints, 1)*(2*window_multiplier*sigma)-window_multiplier*sigma).';
                            Deltas = DiffMat(ValidXfits(k),:);
                            Deltas(Deltas > window_multiplier*sigma | Deltas < -window_multiplier*sigma) = NaN;
                            ind = NaN(1, NumPoints);
                            for l = 1:NumPoints
                                ind(l) = find(min(abs(Deltas - r1(l))) == abs(Deltas - r1(l)));
                            end
                            SmoothedProfiles(ValidXfits(k),:,rep) =  GaussianWeights(ValidXfits(k),ind)*ys(ind, :)/(sum(GaussianWeights(ValidXfits(k),ind)));
                        end
                    end
                    
                    
                    
                    MeanProf = mean(SmoothedProfiles, 3, 'omitnan');
                    SEProf = std(SmoothedProfiles,0,  3, 'omitnan')/sqrt(NumBootstrappedFits);
                    SEProf(MeanProf == 0) = NaN;
                    MeanProf(MeanProf == 0) = NaN;
                    
                    CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).ControlSet.mean(:,:,ch_index, master_index) = MeanProf;
                    CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).ControlSet.se(:,:,ch_index, master_index) = SEProf;
                    CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).ControlSet.counts(:, master_index) = counts;
                    CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).ControlSet.counts_above(:, master_index) = counts_above;
                    CompiledEmbryos.BootstrappedUnivScaledProfiles.(FitTypeStrings{i}).(SetStrings{j}).ControlSet.counts_below(:, master_index) = counts_below;
                    
                end
            end
        end
    end
end




%%

CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');