function CompiledEmbryos =  AddBootstrappedFitProfiles(CompiledEmbryos, exp_index)
AllSetInfo = GetFixedSetPrefixInfo;

AllSetsCombinedEmbryosPath = 'S:/Gabriella/Dropbox/ProteinProfiles/CompiledEmbryos/';

SetLabel = AllSetInfo.SetLabels{exp_index};
OutEmbryoPath = [AllSetsCombinedEmbryosPath, SetLabel];
NChannels = 5;
%%
CompiledEmbryos.UnivBootstrappedScaledProfiles = {};
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP = {};
disp('Need to switch to run over all sets')
for set_index = 1:3
    SetString = ['Set', num2str(set_index)];
    CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString) = {};
    CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).AllProfiles = {};
    CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).AllProfiles.x = ...
        CompiledEmbryos.DubuisEmbryoTimes;
    CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Test = {};
    CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Control = {};
end


%CompiledEmbryos.UnivBootstrappedScaledProfiles.ZeroCorrectedSlideRescaledAvgAP = {};
for i = 1:size(CompiledEmbryos.BootstrappedScaleFactors,1)
    SetString = ['Set', num2str(i)];
    CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).AllProfiles.mean = NaN(size(CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles));
    %CompiledEmbryos.UnivBootstrappedScaledProfiles.ZeroCorrectedSlideRescaledAvgAP.(SetString).mean = NaN(size(CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles));
    for ch_index = 2:NChannels
        CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).AllProfiles.mean(:,:,ch_index) = ...
            CompiledEmbryos.BootstrappedScaleFactors(i, ch_index)*CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + CompiledEmbryos.BootstrappedScaleIntercepts(i, ch_index);
%         CompiledEmbryos.UnivBootstrappedScaledProfiles.ZeroCorrectedSlideRescaledAvgAP.(SetString).mean (:,:,ch_index) = ...
%             CompiledEmbryos.BootstrappedScaleFactors(i, ch_index)*CompiledEmbryos.ZeroCorrectedSlideRescaledDorsalAvgAPProfiles(:,:,ch_index) + CompiledEmbryos.BootstrappedScaleIntercepts(i, ch_index);
    end


end
%%
xfits = 0:.1:70;
Nxfits = length(xfits);
min_2sigma_points = 5;
sigma = 5;
window_multiplier = 3;
NumBootstrappedFits = 1000;
NumPoints = 100;
NumAPbins = size(CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles, 2);
NChannels =size(CompiledEmbryos.FitSlideRescaledDorsalAvgAPProfiles, 3);



%% Bootstrapping the FitSlideRescaledDorsalAvgAPProfiles for the Test Set 
for set_index = 1:3
    SetString = ['Set', num2str(set_index)];
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Test.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Test.se = NaN(Nxfits, NumAPbins, NChannels);
for ch_index = 3
UseTF = CompiledEmbryos.TestSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTF);
ys = CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).AllProfiles.mean(UseTF, :,ch_index);
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
    
    if counts_above(x_index) < 5
        GaussianWeights(x_index,:) = 0;
    end
    
    if counts_below(x_index) < 5
        GaussianWeights(x_index,:) = 0;
    end

end


ValidCountTFs = (counts >= min_2sigma_points & counts_above >= 3 & counts_below >= 3 & xfits >= min(x_sample)+sigma/2 & xfits <= max(x_sample) - sigma/2);
x_indices = 1:Nxfits;
ValidXfits = x_indices(ValidCountTFs);
NValidXfits = length(ValidXfits);

for j = 1:NValidXfits
    %WindowWidthLimits = [xfits(ValidXfits(j))-( min(x_sample)), max(x_sample)- xfits(ValidXfits(j))];
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
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Test.x = xfits;
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Test.mean (:,:,ch_index) = MeanSmoothedProfiles;
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Test.se(:,:,ch_index)  = SmoothedProfileSEs;
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Test.counts = counts;
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Test.counts_above = counts_above;
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Test.counts_abelow = counts_above;

end


%% Bootstrapping the FitSlideRescaledDorsalAvgAPProfiles for the Control Set 
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Control.mean = NaN(Nxfits, NumAPbins, NChannels);
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Control.se = NaN(Nxfits, NumAPbins, NChannels);
for ch_index = 3
UseTF = CompiledEmbryos.ControlSetEmbryos & CompiledEmbryos.IsNC14 & ~isnan(CompiledEmbryos.DubuisEmbryoTimes);
x_sample = CompiledEmbryos.DubuisEmbryoTimes(UseTF);
ys =  CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).AllProfiles.mean(UseTF, :,ch_index);
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
    
    if counts_above(x_index) < 5
        GaussianWeights(x_index,:) = 0;
    end
    
    if counts_below(x_index) < 5
        GaussianWeights(x_index,:) = 0;
    end

end


ValidCountTFs = (counts >= min_2sigma_points & counts_above >= 3 & counts_below >= 3 & xfits >= min(x_sample)+sigma/2 & xfits <= max(x_sample) - sigma/2);
x_indices = 1:Nxfits;
ValidXfits = x_indices(ValidCountTFs);
NValidXfits = length(ValidXfits);

for j = 1:NValidXfits
    %WindowWidthLimits = [xfits(ValidXfits(j))-( min(x_sample)), max(x_sample)- xfits(ValidXfits(j))];
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
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Control.x = xfits;
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Control.mean(:,:,ch_index) = MeanSmoothedProfiles;
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Control.se(:,:,ch_index)  = SmoothedProfileSEs;
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Control.counts = counts;
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Control.counts_above = counts_above;
CompiledEmbryos.UnivBootstrappedScaledProfiles.FitSlideRescaledAvgAP.(SetString).Control.counts_abelow = counts_above;

end
%

%%
end

CEoutpath = [OutEmbryoPath, filesep, 'CompiledEmbryos.mat'];
save(CEoutpath, 'CompiledEmbryos');