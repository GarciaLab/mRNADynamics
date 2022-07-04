function CompiledEmbryos =  AddBootstrappedProfiles(CompiledEmbryos, exp_index)
AllSetInfo = GetFixedSetPrefixInfo;

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];

%%
if ~isfield(CompiledEmbryos, 'BootstrappedProfiles')
    CompiledEmbryos.BootstrappedProfiles = {};
end
if ~isfield(CompiledEmbryos.BootstrappedProfiles, 'FitSlideRescaledDorsalAvgAP')
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP = {};
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.x = {};
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Test = {};
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Control = {};
end

if ~isfield(CompiledEmbryos.BootstrappedProfiles, 'ZeroCorrectedSlideRescaledDorsalAvgAP')
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP = {};
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.x = {};
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Test = {};
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Control = {};
end

if ~isfield(CompiledEmbryos.BootstrappedProfiles, 'SlideRescaledDorsalAvgAP')
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP = {};
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.x = {};
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Test = {};
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Control = {};
end

xfits = 0:.1:70;
Nxfits = length(xfits);
min_2sigma_points = 10;
counts_above_limit = 3;
counts_below_limit = 3;
sigma = 5;
window_multiplier = 3;
NumBootstrappedFits = 500;
NumPoints = 100;
NumAPbins = size(CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles, 2);
NChannels =size(CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles, 3);



%% Bootstrapping the FitSlideRescaledDorsalAvgAPProfiles for the Test Set
if ~isfield(CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Test, 'mean')
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Test.mean = NaN(Nxfits, NumAPbins, NChannels);
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Test.se = NaN(Nxfits, NumAPbins, NChannels);
end
for ch_index = 3
    UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
    x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTF);
    ys = CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(UseTF, :,ch_index);
    BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
    x_sample = x_sample(~BadTF);
    ys = ys(~BadTF,:);
    
    
    
    [x_sample, sort_order] = sort(x_sample);
    ys = ys(sort_order,:);
    
    SmoothedProfiles = NaN(Nxfits, NumAPbins, NumBootstrappedFits);
    
    DiffMat = xfits.'-x_sample;
    GaussianWeights = GetGaussianWeightMat(xfits, x_sample, sigma, window_multiplier*sigma);
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
    
    for j = 1:NValidXfits
        WindowWidthLimits = [xfits(ValidXfits(j))-( min(x_sample)), max(x_sample)- xfits(ValidXfits(j))];
        for rep = 1:NumBootstrappedFits*2
            r1  = (rand(NumPoints, 1)*(2*window_multiplier*sigma)-window_multiplier*sigma).';
            Deltas = DiffMat(ValidXfits(j),:);
            Deltas(Deltas > window_multiplier*sigma | Deltas < -window_multiplier*sigma) = NaN;
            ind = NaN(1, NumPoints);
            for k = 1:NumPoints
                ind(k) = find(min(abs(Deltas - r1(k))) == abs(Deltas - r1(k)));
            end
            SmoothedProfiles(ValidXfits(j),:,rep) =  GaussianWeights(ValidXfits(j),ind)*ys(ind, :)/(sum(GaussianWeights(ValidXfits(j),ind)));
        end
    end
    
    
    MeanSmoothedProfiles = mean(SmoothedProfiles, 3, 'omitnan');
    SmoothedProfileSEs = std(SmoothedProfiles,0,  3, 'omitnan')/sqrt(NumBootstrappedFits);
    
    SmoothedProfileSEs(MeanSmoothedProfiles  == 0) = NaN;
    MeanSmoothedProfiles(MeanSmoothedProfiles  == 0) = NaN;
    %close all
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.x = xfits;
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Test.mean(:,:,ch_index) = MeanSmoothedProfiles;
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Test.se(:,:,ch_index)  = SmoothedProfileSEs;
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Test.counts = counts;
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Test.counts_above = counts_above;
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Test.counts_below = counts_above;
    
end


%% Bootstrapping the FitSlideRescaledDorsalAvgAPProfiles for the Control Set
if ~isfield(CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Control, 'mean')
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Control.mean = NaN(Nxfits, NumAPbins, NChannels);
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Control.se = NaN(Nxfits, NumAPbins, NChannels);
end
for ch_index = 3
    UseTF = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
    x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTF);
    ys = CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(UseTF, :,ch_index);
    BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
    x_sample = x_sample(~BadTF);
    ys = ys(~BadTF,:);
    
    
    
    [x_sample, sort_order] = sort(x_sample);
    ys = ys(sort_order,:);
    
    SmoothedProfiles = NaN(Nxfits, NumAPbins, NumBootstrappedFits);
    
    DiffMat = xfits.'-x_sample;
    GaussianWeights = GetGaussianWeightMat(xfits, x_sample, sigma, window_multiplier*sigma);
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
    
    for j = 1:NValidXfits
        WindowWidthLimits = [xfits(ValidXfits(j))-( min(x_sample)), max(x_sample)- xfits(ValidXfits(j))];
        for rep = 1:NumBootstrappedFits*2
            r1  = (rand(NumPoints, 1)*(2*window_multiplier*sigma)-window_multiplier*sigma).';
            Deltas = DiffMat(ValidXfits(j),:);
            Deltas(Deltas > window_multiplier*sigma | Deltas < -window_multiplier*sigma) = NaN;
            ind = NaN(1, NumPoints);
            for k = 1:NumPoints
                ind(k) = find(min(abs(Deltas - r1(k))) == abs(Deltas - r1(k)));
            end
            SmoothedProfiles(ValidXfits(j),:,rep) =  GaussianWeights(ValidXfits(j),ind)*ys(ind, :)/(sum(GaussianWeights(ValidXfits(j),ind)));
        end
    end
    
    
    
    MeanSmoothedProfiles = mean(SmoothedProfiles, 3, 'omitnan');
    SmoothedProfileSEs = std(SmoothedProfiles,0,  3, 'omitnan')/sqrt(NumBootstrappedFits);
    
    SmoothedProfileSEs(MeanSmoothedProfiles == 0) = NaN;
    MeanSmoothedProfiles(MeanSmoothedProfiles == 0) = NaN;
    %close all
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.x = xfits;
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Control.mean(:,:,ch_index) = MeanSmoothedProfiles;
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Control.se(:,:,ch_index)  = SmoothedProfileSEs;
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Control.counts = counts;
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Control.counts_above = counts_above;
    CompiledEmbryos.BootstrappedProfiles.FitSlideRescaledDorsalAvgAP.Control.counts_below = counts_above;
    
end
%
%% Bootstrapping the ZeroCorrectedSlideRescaledDorsalAvgAP for the Test Set
if ~isfield(CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Test, 'mean')
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Test.mean = NaN(Nxfits, NumAPbins, NChannels);
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Test.se = NaN(Nxfits, NumAPbins, NChannels);
end
for ch_index = [3 5]
    UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
    x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTF);
    ys = CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(UseTF, :,ch_index);
    BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
    x_sample = x_sample(~BadTF);
    ys = ys(~BadTF,:);
    
    
    
    [x_sample, sort_order] = sort(x_sample);
    ys = ys(sort_order,:);
    
    SmoothedProfiles = NaN(Nxfits, NumAPbins, NumBootstrappedFits);
    
    DiffMat = xfits.'-x_sample;
    GaussianWeights = GetGaussianWeightMat(xfits, x_sample, sigma, window_multiplier*sigma);
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
    
    for j = 1:NValidXfits
        WindowWidthLimits = [xfits(ValidXfits(j))-( min(x_sample)), max(x_sample)- xfits(ValidXfits(j))];
        for rep = 1:NumBootstrappedFits*2
            r1  = (rand(NumPoints, 1)*(2*window_multiplier*sigma)-window_multiplier*sigma).';
            Deltas = DiffMat(ValidXfits(j),:);
            Deltas(Deltas > window_multiplier*sigma | Deltas < -window_multiplier*sigma) = NaN;
            ind = NaN(1, NumPoints);
            for k = 1:NumPoints
                ind(k) = find(min(abs(Deltas - r1(k))) == abs(Deltas - r1(k)));
            end
            SmoothedProfiles(ValidXfits(j),:,rep) =  GaussianWeights(ValidXfits(j),ind)*ys(ind, :)/(sum(GaussianWeights(ValidXfits(j),ind)));
        end
    end
    
    
    
    MeanSmoothedProfiles = mean(SmoothedProfiles, 3, 'omitnan');
    SmoothedProfileSEs = std(SmoothedProfiles,0,  3, 'omitnan')/sqrt(NumBootstrappedFits);
    
    SmoothedProfileSEs(MeanSmoothedProfiles== 0) = NaN;
    MeanSmoothedProfiles(MeanSmoothedProfiles  == 0) = NaN;
    %close all
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.x = xfits;
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Test.mean(:,:,ch_index) = MeanSmoothedProfiles;
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Test.se(:,:,ch_index)  = SmoothedProfileSEs;
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Test.counts = counts;
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Test.counts_above = counts_above;
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Test.counts_below = counts_above;
    
end


%% Bootstrapping the FitSlideRescaledDorsalAvgAPProfiles for the Control Set
if ~isfield(CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Control, 'mean')
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Control.mean = NaN(Nxfits, NumAPbins, NChannels);
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Control.se = NaN(Nxfits, NumAPbins, NChannels);
end
for ch_index = [3 5]
    UseTF = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
    x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTF);
    ys = CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(UseTF, :,ch_index);
    BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
    x_sample = x_sample(~BadTF);
    ys = ys(~BadTF,:);
    
    
    
    [x_sample, sort_order] = sort(x_sample);
    ys = ys(sort_order,:);
    
    SmoothedProfiles = NaN(Nxfits, NumAPbins, NumBootstrappedFits);
    
    DiffMat = xfits.'-x_sample;
    GaussianWeights = GetGaussianWeightMat(xfits, x_sample, sigma, window_multiplier*sigma);
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
    
    for j = 1:NValidXfits
        WindowWidthLimits = [xfits(ValidXfits(j))-( min(x_sample)), max(x_sample)- xfits(ValidXfits(j))];
        for rep = 1:NumBootstrappedFits*2
            r1  = (rand(NumPoints, 1)*(2*window_multiplier*sigma)-window_multiplier*sigma).';
            Deltas = DiffMat(ValidXfits(j),:);
            Deltas(Deltas > window_multiplier*sigma | Deltas < -window_multiplier*sigma) = NaN;
            ind = NaN(1, NumPoints);
            for k = 1:NumPoints
                ind(k) = find(min(abs(Deltas - r1(k))) == abs(Deltas - r1(k)));
            end
            SmoothedProfiles(ValidXfits(j),:,rep) =  GaussianWeights(ValidXfits(j),ind)*ys(ind, :)/(sum(GaussianWeights(ValidXfits(j),ind)));
        end
    end
    
    
    MeanSmoothedProfiles = mean(SmoothedProfiles, 3, 'omitnan');
    SmoothedProfileSEs = std(SmoothedProfiles,0,  3, 'omitnan')/sqrt(NumBootstrappedFits);
    
    SmoothedProfileSEs(MeanSmoothedProfiles == 0) = NaN;
    MeanSmoothedProfiles(MeanSmoothedProfiles  == 0) = NaN;
    %close all
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.x = xfits;
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Control.mean(:,:,ch_index) = MeanSmoothedProfiles;
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Control.se(:,:,ch_index)  = SmoothedProfileSEs;
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Control.counts = counts;
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Control.counts_above = counts_above;
    CompiledEmbryos.BootstrappedProfiles.ZeroCorrectedSlideRescaledDorsalAvgAP.Control.counts_below = counts_above;
end

%% Bootstrapping the SlideRescaledDorsalAvgAP for the Test Set
if ~isfield(CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Test, 'mean')
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Test.mean = NaN(Nxfits, NumAPbins, NChannels);
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Test.se = NaN(Nxfits, NumAPbins, NChannels);
end
for ch_index = [3 5]
    UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
    x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTF);
    ys = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(UseTF, :,ch_index);
    BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
    x_sample = x_sample(~BadTF);
    ys = ys(~BadTF,:);
    
    
    
    [x_sample, sort_order] = sort(x_sample);
    ys = ys(sort_order,:);
    
    SmoothedProfiles = NaN(Nxfits, NumAPbins, NumBootstrappedFits);
    
    DiffMat = xfits.'-x_sample;
    GaussianWeights = GetGaussianWeightMat(xfits, x_sample, sigma, window_multiplier*sigma);
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
    
    for j = 1:NValidXfits
        WindowWidthLimits = [xfits(ValidXfits(j))-( min(x_sample)), max(x_sample)- xfits(ValidXfits(j))];
        for rep = 1:NumBootstrappedFits*2
            r1  = (rand(NumPoints, 1)*(2*window_multiplier*sigma)-window_multiplier*sigma).';
            Deltas = DiffMat(ValidXfits(j),:);
            Deltas(Deltas > window_multiplier*sigma | Deltas < -window_multiplier*sigma) = NaN;
            ind = NaN(1, NumPoints);
            for k = 1:NumPoints
                ind(k) = find(min(abs(Deltas - r1(k))) == abs(Deltas - r1(k)));
            end
            SmoothedProfiles(ValidXfits(j),:,rep) =  GaussianWeights(ValidXfits(j),ind)*ys(ind, :)/(sum(GaussianWeights(ValidXfits(j),ind)));
        end
    end
    
    
    MeanSmoothedProfiles = mean(SmoothedProfiles, 3, 'omitnan');
    SmoothedProfileSEs = std(SmoothedProfiles,0,  3, 'omitnan')/sqrt(NumBootstrappedFits);
    
    SmoothedProfileSEs(MeanSmoothedProfiles== 0) = NaN;
    MeanSmoothedProfiles(MeanSmoothedProfiles  == 0) = NaN;
    %close all
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.x = xfits;
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Test.mean(:,:,ch_index) = MeanSmoothedProfiles;
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Test.se(:,:,ch_index)  = SmoothedProfileSEs;
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Test.counts = counts;
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Test.counts_above = counts_above;
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Test.counts_below = counts_above;
    
end


%% Bootstrapping the SlideRescaledDorsalAvgAP for the Control Set
if ~isfield(CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Control, 'mean')
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Control.mean = NaN(Nxfits, NumAPbins, NChannels);
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Control.se = NaN(Nxfits, NumAPbins, NChannels);
end
for ch_index = [3 5]
    UseTF = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
    x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTF);
    ys = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles(UseTF, :,ch_index);
    BadTF =  sum(isnan(ys), 2).' >min( sum(isnan(ys), 2).');
    x_sample = x_sample(~BadTF);
    ys = ys(~BadTF,:);
    
    
    
    [x_sample, sort_order] = sort(x_sample);
    ys = ys(sort_order,:);
    
    SmoothedProfiles = NaN(Nxfits, NumAPbins, NumBootstrappedFits);
    
    DiffMat = xfits.'-x_sample;
    GaussianWeights = GetGaussianWeightMat(xfits, x_sample, sigma, window_multiplier*sigma);
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
    
    for j = 1:NValidXfits
        WindowWidthLimits = [xfits(ValidXfits(j))-( min(x_sample)), max(x_sample)- xfits(ValidXfits(j))];
        for rep = 1:NumBootstrappedFits*2
            r1  = (rand(NumPoints, 1)*(2*window_multiplier*sigma)-window_multiplier*sigma).';
            Deltas = DiffMat(ValidXfits(j),:);
            Deltas(Deltas > window_multiplier*sigma | Deltas < -window_multiplier*sigma) = NaN;
            ind = NaN(1, NumPoints);
            for k = 1:NumPoints
                ind(k) = find(min(abs(Deltas - r1(k))) == abs(Deltas - r1(k)));
            end
            SmoothedProfiles(ValidXfits(j),:,rep) =  GaussianWeights(ValidXfits(j),ind)*ys(ind, :)/(sum(GaussianWeights(ValidXfits(j),ind)));
        end
    end
    
    
    
    MeanSmoothedProfiles = mean(SmoothedProfiles, 3, 'omitnan');
    SmoothedProfileSEs = std(SmoothedProfiles,0,  3, 'omitnan')/sqrt(NumBootstrappedFits);
    
    SmoothedProfileSEs(MeanSmoothedProfiles == 0) = NaN;
    MeanSmoothedProfiles(MeanSmoothedProfiles  == 0) = NaN;
    %close all
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.x = xfits;
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Control.mean(:,:,ch_index) = MeanSmoothedProfiles;
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Control.se(:,:,ch_index)  = SmoothedProfileSEs;
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Control.counts = counts;
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Control.counts_above = counts_above;
    CompiledEmbryos.BootstrappedProfiles.SlideRescaledDorsalAvgAP.Control.counts_below = counts_above;
end
%%


CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');