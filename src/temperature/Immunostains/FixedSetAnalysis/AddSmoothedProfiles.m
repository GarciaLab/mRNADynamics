function CompiledEmbryos = AddSmoothedProfiles(CompiledEmbryos, sigma_deltafc, sigma_time_multiplier, width_sig_multiplier, min_points_singleside, min_2sigma_points)%, exp_index)
%%
if ~exist('sigma_deltafc', 'var')
    sigma_deltafc = 2.5;
end
if ~exist('sigma_time_multiplier', 'var')
    sigma_time_multiplier = 1/10;
end
    
if ~exist('width_sig_multiplier', 'var')
    width_sig_multiplier = 3;
end

if ~exist('min_points_singleside', 'var')
    min_points_singleside = 1;
end

if ~exist('min_2sigma_points', 'var')
    min_2sigma_points = 3;
end

if ~exist('sigma_time', 'var')
    sigma_time = 4;
end

APbins = 0:0.025:1;
NumAPbins = length(APbins);
NarrowAPbins = 0:0.0125:1;
NumNarrowAPbins = length(NarrowAPbins);

CompiledEmbryos.AllDeltaValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan([CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].');% &...
        %CompiledEmbryos.FixCorrectedDeltaFC_um(:,1).' > 3;
DorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles;
NChannels = size(DorsalProfiles, 3);
y = CompiledEmbryos.FixCorrectedDeltaFC_um.mean(CompiledEmbryos.AllDeltaValidProfilesTestTF);
% ThreeMinDeltas = mink(y, 3);
% minx = floor(ThreeMinDeltas(end)/.1)*.1;
x = 2:1:45;
sigma = sigma_deltafc;
DiffMat = x.'-y;
[nn_idx, nn_D] = knnsearch(y.', x.', 'K', 50);



GaussianWeights = GetGaussianWeightMat(x, y, sigma, (width_sig_multiplier)*sigma);
keepx = ones(1, length(x), 'logical');
counts_within_halfsigma = zeros(1, length(x));
counts_within_1sigma = zeros(1, length(x));
counts_within_2sigma = zeros(1, length(x));
counts_within_3sigma = zeros(1, length(x));
counts_within_4sigma = zeros(1, length(x));
for x_index = 1:length(x)
    Num2SigmaPoints = sum(nn_D(x_index,:) <= 2*sigma_deltafc);
    
    if Num2SigmaPoints < min_2sigma_points
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    if sum(DiffMat(x_index, nn_idx(x_index,:)) > 0) < min_points_singleside | sum(DiffMat(x_index, nn_idx(x_index,:)) < 0) < min_points_singleside 
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    
    counts_within_halfsigma(x_index) = sum(nn_D(x_index,:) <= sigma_deltafc/2);
    counts_within_1sigma(x_index) = sum(nn_D(x_index,:) <= sigma_deltafc);
    counts_within_2sigma(x_index) = sum(nn_D(x_index,:) <= 2*sigma_deltafc & nn_D(x_index,:) > sigma_deltafc);
    counts_within_3sigma(x_index) = sum(nn_D(x_index,:) <= 3*sigma_deltafc & nn_D(x_index,:) > 2*sigma_deltafc);
    counts_within_4sigma(x_index) = sum(nn_D(x_index,:) <= 4*sigma_deltafc & nn_D(x_index,:) > 3*sigma_deltafc);
    
    GaussianWeights(x_index, ~ismember(1:length(y), nn_idx(x_index,:))) = 0;
   
end
CompiledEmbryos.FixCorrectedSmoothedProfiles = {};
CompiledEmbryos.FixCorrectedSmoothedProfiles.TestCounts = {};
CompiledEmbryos.FixCorrectedSmoothedProfiles.TestCounts.counts_within_1sigma = counts_within_1sigma;
CompiledEmbryos.FixCorrectedSmoothedProfiles.TestCounts.counts_within_2sigma = counts_within_2sigma;
CompiledEmbryos.FixCorrectedSmoothedProfiles.TestCounts.counts_within_3sigma = counts_within_3sigma;
CompiledEmbryos.FixCorrectedSmoothedProfiles.TestCounts.counts_within_4sigma = counts_within_4sigma;

CompiledEmbryos.FixCorrectedSmoothedDeltaFCs = x;
CompiledEmbryos.FixCorrectedSmoothedAvgAPProfiles = {};
CompiledEmbryos.FixCorrectedSmoothedAvgAPProfiles.Test = NaN(length(x), NumAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:length(APbins)
        BinData = DorsalProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.FixCorrectedSmoothedAvgAPProfiles.Test(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end

CompiledEmbryos.FixCorrectedSmoothedAvgNarrowAPProfiles ={};
CompiledEmbryos.FixCorrectedSmoothedAvgNarrowAPProfiles.Test = NaN(length(x), NumNarrowAPbins,  NChannels);
NarrowDorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles;
for ch_index = 2:NChannels % 2:5
    for RefBin = 1:NumNarrowAPbins
        BinData = NarrowDorsalProfiles(CompiledEmbryos.AllDeltaValidProfilesTestTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.FixCorrectedSmoothedAvgNarrowAPProfiles.Test(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end


CompiledEmbryos.AllDeltaValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan([CompiledEmbryos.FixCorrectedDeltaFC_um.mean(:)].');% &...
    
y = CompiledEmbryos.FixCorrectedDeltaFC_um.mean(CompiledEmbryos.AllDeltaValidProfilesControlTF);
% ThreeMinDeltas = mink(y, 3);
% minx = floor(ThreeMinDeltas(end)/.1)*.1;
x = 2:1:45;
sigma = sigma_deltafc;



DiffMat = x.'-y;
[nn_idx, nn_D] = knnsearch(y.', x.', 'K', 50);
GaussianWeights = GetGaussianWeightMat(x, y, sigma, (width_sig_multiplier)*sigma);
keepx = ones(1, length(x), 'logical');
counts_within_1sigma = zeros(1, length(x));
counts_within_2sigma = zeros(1, length(x));
counts_within_3sigma = zeros(1, length(x));
counts_within_4sigma = zeros(1, length(x));
for x_index = 1:length(x)
    Num2SigmaPoints = sum(nn_D(x_index,:) <= 2*sigma_deltafc);
    
    if Num2SigmaPoints < min_2sigma_points
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    if sum(DiffMat(x_index, nn_idx(x_index,:)) > 0) < min_points_singleside | sum(DiffMat(x_index, nn_idx(x_index,:)) < 0) < min_points_singleside 
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    
    counts_within_1sigma(x_index) = sum(nn_D(x_index,:) <= sigma_deltafc);
    counts_within_2sigma(x_index) = sum(nn_D(x_index,:) <= 2*sigma_deltafc & nn_D(x_index,:) > sigma_deltafc);
    counts_within_3sigma(x_index) = sum(nn_D(x_index,:) <= 3*sigma_deltafc & nn_D(x_index,:) > 2*sigma_deltafc);
    counts_within_4sigma(x_index) = sum(nn_D(x_index,:) <= 4*sigma_deltafc & nn_D(x_index,:) > 3*sigma_deltafc);
    
    GaussianWeights(x_index, ~ismember(1:length(y), nn_idx(x_index,:))) = 0;
   
end
CompiledEmbryos.FixCorrectedSmoothedProfiles.ControlCounts = {};
CompiledEmbryos.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_1sigma = counts_within_1sigma;
CompiledEmbryos.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_2sigma = counts_within_2sigma;
CompiledEmbryos.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_3sigma = counts_within_3sigma;
CompiledEmbryos.FixCorrectedSmoothedProfiles.ControlCounts.counts_within_4sigma = counts_within_4sigma;

CompiledEmbryos.FixCorrectedSmoothedAvgAPProfiles.Control = NaN(length(x), NumAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:length(APbins)
        BinData = DorsalProfiles(CompiledEmbryos.AllDeltaValidProfilesControlTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.FixCorrectedSmoothedAvgAPProfiles.Control(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end

CompiledEmbryos.FixCorrectedSmoothedAvgNarrowAPProfiles.Control = NaN(length(x), NumNarrowAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:NumNarrowAPbins
        BinData = NarrowDorsalProfiles(CompiledEmbryos.AllDeltaValidProfilesControlTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.FixCorrectedSmoothedAvgNarrowAPProfiles.Control(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end


        %CompiledEmbryos.FixCorrectedDeltaFC_um(:,1).' > 3;

%%

CompiledEmbryos.AllDubuisValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan([CompiledEmbryos.DubuisEmbryoTimes]);% &...
y = CompiledEmbryos.DubuisEmbryoTimes(CompiledEmbryos.AllDubuisValidProfilesTestTF );
% ThreeMinDeltas = mink(y, 3);
% minx = floor(ThreeMinDeltas(end)/.1)*.1;
MaxT = ceil(max([CompiledEmbryos.DubuisEmbryoTimes(:)]));
MaxT = 70;
x = 0:1:MaxT;
DorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgAPProfiles;

sigma = sigma_time;%floor(MaxT*sigma_time_multiplier);

DiffMat = x.'-y;
[nn_idx, nn_D] = knnsearch(y.', x.', 'K', 50);
GaussianWeights = GetGaussianWeightMat(x, y, sigma, (width_sig_multiplier)*sigma);
keepx = ones(1, length(x), 'logical');
counts_within_halfsigma = zeros(1, length(x));
counts_within_1sigma = zeros(1, length(x));
counts_within_2sigma = zeros(1, length(x));
counts_within_3sigma = zeros(1, length(x));
counts_within_4sigma = zeros(1, length(x));
for x_index = 1:length(x)
    Num2SigmaPoints = sum(nn_D(x_index,:) <= 2*sigma);
    
    if Num2SigmaPoints < min_2sigma_points
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    if sum(DiffMat(x_index, nn_idx(x_index,:)) > 0) < min_points_singleside | sum(DiffMat(x_index, nn_idx(x_index,:)) < 0) < min_points_singleside 
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    
    counts_within_halfsigma(x_index)= sum(nn_D(x_index,:) <= sigma/2);
    counts_within_1sigma(x_index) = sum(nn_D(x_index,:) <= sigma);
    counts_within_2sigma(x_index) = sum(nn_D(x_index,:) <= 2*sigma & nn_D(x_index,:) > sigma);
    counts_within_3sigma(x_index) = sum(nn_D(x_index,:) <= 3*sigma & nn_D(x_index,:) > 2*sigma);
    counts_within_4sigma(x_index) = sum(nn_D(x_index,:) <= 4*sigma & nn_D(x_index,:) > 3*sigma);
    
    GaussianWeights(x_index, ~ismember(1:length(y), nn_idx(x_index,:))) = 0;
   
end
CompiledEmbryos.DubuisTimesSmoothedProfiles = {};
CompiledEmbryos.DubuisTimesSmoothedProfiles.TestCounts = {};
CompiledEmbryos.DubuisTimesSmoothedProfiles.TestCounts.counts_within_halfsigma = counts_within_halfsigma;
CompiledEmbryos.DubuisTimesSmoothedProfiles.TestCounts.counts_within_1sigma = counts_within_1sigma;
CompiledEmbryos.DubuisTimesSmoothedProfiles.TestCounts.counts_within_2sigma = counts_within_2sigma;
CompiledEmbryos.DubuisTimesSmoothedProfiles.TestCounts.counts_within_3sigma = counts_within_3sigma;
CompiledEmbryos.DubuisTimesSmoothedProfiles.TestCounts.counts_within_4sigma = counts_within_4sigma;

CompiledEmbryos.DubuisSmoothedTimes = x;
CompiledEmbryos.DubuisTimeSmoothedAvgAPProfiles = {};
CompiledEmbryos.DubuisTimeSmoothedAvgAPProfiles.Test = NaN(length(x), NumAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:length(APbins)
        BinData = DorsalProfiles(CompiledEmbryos.AllDubuisValidProfilesTestTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.DubuisTimeSmoothedAvgAPProfiles.Test(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end


NarrowDorsalProfiles = CompiledEmbryos.SlideRescaledDorsalAvgNarrowAPProfiles;
CompiledEmbryos.DubuisTimeSmoothedAvgNarrowAPProfiles = {};
CompiledEmbryos.DubuisTimeSmoothedAvgNarrowAPProfiles.Test = NaN(length(x), NumNarrowAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:NumNarrowAPbins
        BinData = NarrowDorsalProfiles(CompiledEmbryos.AllDubuisValidProfilesTestTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.DubuisTimeSmoothedAvgNarrowAPProfiles.Test(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end

CompiledEmbryos.AllDubuisValidProfilesControlTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan([CompiledEmbryos.DubuisEmbryoTimes]);% &...
y = CompiledEmbryos.DubuisEmbryoTimes(CompiledEmbryos.AllDubuisValidProfilesControlTF );
% ThreeMinDeltas = mink(y, 3);
% minx = floor(ThreeMinDeltas(end)/.1)*.1;
MaxT = ceil(max([CompiledEmbryos.DubuisEmbryoTimes(:)]));
MaxT = 70;
x = 0:1:MaxT;

DiffMat = x.'-y;
[nn_idx, nn_D] = knnsearch(y.', x.', 'K', 50);
sigma = sigma_time;%floor(MaxT*sigma_time_multiplier);
GaussianWeights = GetGaussianWeightMat(x, y, sigma, (width_sig_multiplier)*sigma);
keepx = ones(1, length(x), 'logical');
counts_within_halfsigma = zeros(1, length(x));
counts_within_1sigma = zeros(1, length(x));
counts_within_2sigma = zeros(1, length(x));
counts_within_3sigma = zeros(1, length(x));
counts_within_4sigma = zeros(1, length(x));
for x_index = 1:length(x)
    Num2SigmaPoints = sum(nn_D(x_index,:) <= 2*sigma);
    
    if Num2SigmaPoints < min_2sigma_points
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    if sum(DiffMat(x_index, nn_idx(x_index,:)) > 0) < min_points_singleside | sum(DiffMat(x_index, nn_idx(x_index,:)) < 0) < min_points_singleside 
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    
    counts_within_halfsigma(x_index)= sum(nn_D(x_index,:) <= sigma/2);
    counts_within_1sigma(x_index) = sum(nn_D(x_index,:) <= sigma);
    counts_within_2sigma(x_index) = sum(nn_D(x_index,:) <= 2*sigma & nn_D(x_index,:) > sigma);
    counts_within_3sigma(x_index) = sum(nn_D(x_index,:) <= 3*sigma & nn_D(x_index,:) > 2*sigma);
    counts_within_4sigma(x_index) = sum(nn_D(x_index,:) <= 4*sigma & nn_D(x_index,:) > 3*sigma);
    
    GaussianWeights(x_index, ~ismember(1:length(y), nn_idx(x_index,:))) = 0;
   
end

CompiledEmbryos.DubuisTimesSmoothedProfiles.ControlCounts = {};
CompiledEmbryos.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_halfsigma = counts_within_halfsigma;
CompiledEmbryos.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_1sigma = counts_within_1sigma;
CompiledEmbryos.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_2sigma = counts_within_2sigma;
CompiledEmbryos.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_3sigma = counts_within_3sigma;
CompiledEmbryos.DubuisTimesSmoothedProfiles.ControlCounts.counts_within_4sigma = counts_within_4sigma;

CompiledEmbryos.DubuisTimeSmoothedAvgAPProfiles.Control = NaN(length(x), NumAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:length(APbins)
        BinData = DorsalProfiles(CompiledEmbryos.AllDubuisValidProfilesControlTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.DubuisTimeSmoothedAvgAPProfiles.Control(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end

CompiledEmbryos.DubuisTimeSmoothedAvgNarrowAPProfiles.Control = NaN(length(x), NumNarrowAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:NumNarrowAPbins
        BinData = NarrowDorsalProfiles(CompiledEmbryos.AllDubuisValidProfilesControlTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.DubuisTimeSmoothedAvgNarrowAPProfiles.Control(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end

%% Add yw25CsquishedTimes

CompiledEmbryos.Allyw25CValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan([CompiledEmbryos.yw25CEmbryoTimes(:)].');% &...
y = CompiledEmbryos.yw25CEmbryoTimes(CompiledEmbryos.Allyw25CValidProfilesTestTF );
% ThreeMinDeltas = mink(y, 3);
% minx = floor(ThreeMinDeltas(end)/.1)*.1;
MaxT = ceil(max([CompiledEmbryos.yw25CEmbryoTimes(:)]));
MaxT = 70;
x = 0:1:MaxT;
sigma = sigma_time;%floor(MaxT*sigma_time_multiplier);
DiffMat = x.'-y;
[nn_idx, nn_D] = knnsearch(y.', x.', 'K', 25);
GaussianWeights = GetGaussianWeightMat(x, y, sigma, (width_sig_multiplier)*sigma);
keepx = ones(1, length(x), 'logical');
counts_within_1sigma = zeros(1, length(x));
counts_within_2sigma = zeros(1, length(x));
counts_within_3sigma = zeros(1, length(x));
counts_within_4sigma = zeros(1, length(x));
for x_index = 1:length(x)
    Num2SigmaPoints = sum(nn_D(x_index,:) <= 2*sigma_deltafc);
    
    if Num2SigmaPoints < min_2sigma_points
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    if sum(DiffMat(x_index, nn_idx(x_index,:)) > 0) < min_points_singleside | sum(DiffMat(x_index, nn_idx(x_index,:)) < 0) < min_points_singleside 
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    
    counts_within_1sigma(x_index) = sum(nn_D(x_index,:) <= sigma_deltafc);
    counts_within_2sigma(x_index) = sum(nn_D(x_index,:) <= 2*sigma_deltafc & nn_D(x_index,:) > sigma_deltafc);
    counts_within_3sigma(x_index) = sum(nn_D(x_index,:) <= 3*sigma_deltafc & nn_D(x_index,:) > 2*sigma_deltafc);
    counts_within_4sigma(x_index) = sum(nn_D(x_index,:) <= 4*sigma_deltafc & nn_D(x_index,:) > 3*sigma_deltafc);
    
    GaussianWeights(x_index, ~ismember(1:length(y), nn_idx(x_index,:))) = 0;
   
end
CompiledEmbryos.yw25CTimesSmoothedProfiles = {};
CompiledEmbryos.yw25CTimesSmoothedProfiles.TestCounts = {};
CompiledEmbryos.yw25CTimesSmoothedProfiles.TestCounts.counts_within_1sigma = counts_within_1sigma;
CompiledEmbryos.yw25CTimesSmoothedProfiles.TestCounts.counts_within_2sigma = counts_within_2sigma;
CompiledEmbryos.yw25CTimesSmoothedProfiles.TestCounts.counts_within_3sigma = counts_within_3sigma;
CompiledEmbryos.yw25CTimesSmoothedProfiles.TestCounts.counts_within_4sigma = counts_within_4sigma;


CompiledEmbryos.yw25CSmoothedTimes= x;
CompiledEmbryos.yw25CTimeSmoothedAvgAPProfiles = {};
CompiledEmbryos.yw25CTimeSmoothedAvgAPProfiles.Test = NaN(length(x), NumAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:length(APbins)
        BinData = DorsalProfiles(CompiledEmbryos.Allyw25CValidProfilesTestTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.yw25CTimeSmoothedAvgAPProfiles.Test(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end

CompiledEmbryos.yw25CTimeSmoothedAvgNarrowAPProfiles = {};
CompiledEmbryos.yw25CTimeSmoothedAvgNarrowAPProfiles.Test = NaN(length(x), NumNarrowAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:NumNarrowAPbins
        BinData = NarrowDorsalProfiles(CompiledEmbryos.Allyw25CValidProfilesTestTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.yw25CTimeSmoothedAvgNarrowAPProfiles.Test(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end

CompiledEmbryos.Allyw25CValidProfilesControlTF  = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan([CompiledEmbryos.DubuisEmbryoTimes(:)].');% &...
y = CompiledEmbryos.yw25CEmbryoTimes(CompiledEmbryos.Allyw25CValidProfilesControlTF );
% ThreeMinDeltas = mink(y, 3);
% minx = floor(ThreeMinDeltas(end)/.1)*.1;
MaxT = ceil(max([CompiledEmbryos.yw25CEmbryoTimes(:)]));
MaxT = 70;
x = 0:1:MaxT;
sigma = sigma_time;%floor(MaxT*sigma_time_multiplier);
DiffMat = x.'-y;
[nn_idx, nn_D] = knnsearch(y.', x.', 'K', 25);
GaussianWeights = GetGaussianWeightMat(x, y, sigma, (width_sig_multiplier)*sigma);
keepx = ones(1, length(x), 'logical');
counts_within_1sigma = zeros(1, length(x));
counts_within_2sigma = zeros(1, length(x));
counts_within_3sigma = zeros(1, length(x));
counts_within_4sigma = zeros(1, length(x));
for x_index = 1:length(x)
    Num2SigmaPoints = sum(nn_D(x_index,:) <= 2*sigma_deltafc);
    
    if Num2SigmaPoints < min_2sigma_points
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    if sum(DiffMat(x_index, nn_idx(x_index,:)) > 0) < min_points_singleside | sum(DiffMat(x_index, nn_idx(x_index,:)) < 0) < min_points_singleside 
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    
    counts_within_1sigma(x_index) = sum(nn_D(x_index,:) <= sigma_deltafc);
    counts_within_2sigma(x_index) = sum(nn_D(x_index,:) <= 2*sigma_deltafc & nn_D(x_index,:) > sigma_deltafc);
    counts_within_3sigma(x_index) = sum(nn_D(x_index,:) <= 3*sigma_deltafc & nn_D(x_index,:) > 2*sigma_deltafc);
    counts_within_4sigma(x_index) = sum(nn_D(x_index,:) <= 4*sigma_deltafc & nn_D(x_index,:) > 3*sigma_deltafc);
    
    GaussianWeights(x_index, ~ismember(1:length(y), nn_idx(x_index,:))) = 0;
   
end

CompiledEmbryos.yw25CTimesSmoothedProfiles.ControlCounts = {};
CompiledEmbryos.yw25CTimesSmoothedProfiles.ControlCounts.counts_within_1sigma = counts_within_1sigma;
CompiledEmbryos.yw25CTimesSmoothedProfiles.ControlCounts.counts_within_2sigma = counts_within_2sigma;
CompiledEmbryos.yw25CTimesSmoothedProfiles.ControlCounts.counts_within_3sigma = counts_within_3sigma;
CompiledEmbryos.yw25CTimesSmoothedProfiles.ControlCounts.counts_within_4sigma = counts_within_4sigma;


CompiledEmbryos.yw25CTimeSmoothedAvgAPProfiles.Control = NaN(length(x), NumAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:length(APbins)
        BinData = DorsalProfiles(CompiledEmbryos.Allyw25CValidProfilesControlTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.yw25CTimeSmoothedAvgAPProfiles.Control(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end

CompiledEmbryos.yw25CTimeSmoothedAvgNarrowAPProfiles.Control = NaN(length(x), NumNarrowAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:NumNarrowAPbins
        BinData = NarrowDorsalProfiles(CompiledEmbryos.Allyw25CValidProfilesControlTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.yw25CTimeSmoothedAvgNarrowAPProfiles.Control(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end


%% Add HisRFP25CsquishedTimes

CompiledEmbryos.AllHisRFP25CValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan([CompiledEmbryos.HisRFP25CEmbryoTimes(:)].');% &...
y = CompiledEmbryos.HisRFP25CEmbryoTimes(CompiledEmbryos.AllHisRFP25CValidProfilesTestTF );
% ThreeMinDeltas = mink(y, 3);
% minx = floor(ThreeMinDeltas(end)/.1)*.1;
MaxT = ceil(max([CompiledEmbryos.HisRFP25CEmbryoTimes(:)]));
MaxT = 70;
x = 0:1:MaxT;
sigma = sigma_time;%floor(MaxT*sigma_time_multiplier);
DiffMat = x.'-y;
[nn_idx, nn_D] = knnsearch(y.', x.', 'K', 25);
GaussianWeights = GetGaussianWeightMat(x, y, sigma, (width_sig_multiplier)*sigma);
keepx = ones(1, length(x), 'logical');
counts_within_1sigma = zeros(1, length(x));
counts_within_2sigma = zeros(1, length(x));
counts_within_3sigma = zeros(1, length(x));
counts_within_4sigma = zeros(1, length(x));
for x_index = 1:length(x)
    Num2SigmaPoints = sum(nn_D(x_index,:) <= 2*sigma_deltafc);
    
    if Num2SigmaPoints < min_2sigma_points
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    if sum(DiffMat(x_index, nn_idx(x_index,:)) > 0) < min_points_singleside | sum(DiffMat(x_index, nn_idx(x_index,:)) < 0) < min_points_singleside 
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    
    counts_within_1sigma(x_index) = sum(nn_D(x_index,:) <= sigma_deltafc);
    counts_within_2sigma(x_index) = sum(nn_D(x_index,:) <= 2*sigma_deltafc & nn_D(x_index,:) > sigma_deltafc);
    counts_within_3sigma(x_index) = sum(nn_D(x_index,:) <= 3*sigma_deltafc & nn_D(x_index,:) > 2*sigma_deltafc);
    counts_within_4sigma(x_index) = sum(nn_D(x_index,:) <= 4*sigma_deltafc & nn_D(x_index,:) > 3*sigma_deltafc);
    
    GaussianWeights(x_index, ~ismember(1:length(y), nn_idx(x_index,:))) = 0;
   
end
CompiledEmbryos.HisRFPTimesSmoothedProfiles = {};
CompiledEmbryos.HisRFPTimesSmoothedProfiles.TestCounts = {};
CompiledEmbryos.HisRFPTimesSmoothedProfiles.TestCounts.counts_within_1sigma = counts_within_1sigma;
CompiledEmbryos.HisRFPTimesSmoothedProfiles.TestCounts.counts_within_2sigma = counts_within_2sigma;
CompiledEmbryos.HisRFPTimesSmoothedProfiles.TestCounts.counts_within_3sigma = counts_within_3sigma;
CompiledEmbryos.HisRFPTimesSmoothedProfiles.TestCounts.counts_within_4sigma = counts_within_4sigma;


CompiledEmbryos.HisRFP25CSmoothedTimes = x;
CompiledEmbryos.HisRFP25CTimeSmoothedAvgAPProfiles = {};
CompiledEmbryos.HisRFP25CTimeSmoothedAvgAPProfiles.Test = NaN(length(x), NumAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:length(APbins)
        BinData = DorsalProfiles(CompiledEmbryos.AllHisRFP25CValidProfilesTestTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.HisRFP25CTimeSmoothedAvgAPProfiles.Test(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end

CompiledEmbryos.HisRFP25CTimeSmoothedAvgNarrowAPProfiles = {};
CompiledEmbryos.HisRFP25CTimeSmoothedAvgNarrowAPProfiles.Test = NaN(length(x), NumNarrowAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:NumNarrowAPbins
        BinData = NarrowDorsalProfiles(CompiledEmbryos.AllHisRFP25CValidProfilesTestTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.HisRFP25CTimeSmoothedAvgNarrowAPProfiles.Test(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end

CompiledEmbryos.AllHisRFP25CValidProfilesControlTF  = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.ControlSetEmbryos & ~isnan([CompiledEmbryos.HisRFP25CEmbryoTimes(:)].');% &...
y = CompiledEmbryos.HisRFP25CEmbryoTimes(CompiledEmbryos.AllHisRFP25CValidProfilesControlTF );
% ThreeMinDeltas = mink(y, 3);
% minx = floor(ThreeMinDeltas(end)/.1)*.1;
MaxT = ceil(max([CompiledEmbryos.HisRFP25CEmbryoTimes(:)]));
MaxT = 70;
x = 0:1:MaxT;
sigma = sigma_time;%floor(MaxT*sigma_time_multiplier);
DiffMat = x.'-y;
[nn_idx, nn_D] = knnsearch(y.', x.', 'K', 25);
GaussianWeights = GetGaussianWeightMat(x, y, sigma, (width_sig_multiplier)*sigma);
keepx = ones(1, length(x), 'logical');
counts_within_1sigma = zeros(1, length(x));
counts_within_2sigma = zeros(1, length(x));
counts_within_3sigma = zeros(1, length(x));
counts_within_4sigma = zeros(1, length(x));
for x_index = 1:length(x)
    Num2SigmaPoints = sum(nn_D(x_index,:) <= 2*sigma_deltafc);
    
    if Num2SigmaPoints < min_2sigma_points
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    if sum(DiffMat(x_index, nn_idx(x_index,:)) > 0) < min_points_singleside | sum(DiffMat(x_index, nn_idx(x_index,:)) < 0) < min_points_singleside 
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    
    counts_within_1sigma(x_index) = sum(nn_D(x_index,:) <= sigma_deltafc);
    counts_within_2sigma(x_index) = sum(nn_D(x_index,:) <= 2*sigma_deltafc & nn_D(x_index,:) > sigma_deltafc);
    counts_within_3sigma(x_index) = sum(nn_D(x_index,:) <= 3*sigma_deltafc & nn_D(x_index,:) > 2*sigma_deltafc);
    counts_within_4sigma(x_index) = sum(nn_D(x_index,:) <= 4*sigma_deltafc & nn_D(x_index,:) > 3*sigma_deltafc);
    
    GaussianWeights(x_index, ~ismember(1:length(y), nn_idx(x_index,:))) = 0;
   
end

CompiledEmbryos.HisRFPTimesSmoothedProfiles.ControlCounts = {};
CompiledEmbryos.HisRFPTimesSmoothedProfiles.ControlCounts.counts_within_1sigma = counts_within_1sigma;
CompiledEmbryos.HisRFPTimesSmoothedProfiles.ControlCounts.counts_within_2sigma = counts_within_2sigma;
CompiledEmbryos.HisRFPTimesSmoothedProfiles.ControlCounts.counts_within_3sigma = counts_within_3sigma;
CompiledEmbryos.HisRFPTimesSmoothedProfiles.ControlCounts.counts_within_4sigma = counts_within_4sigma;

CompiledEmbryos.HisRFP25CTimeSmoothedAvgAPProfiles.Control = NaN(length(x), NumAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:length(APbins)
        BinData = DorsalProfiles(CompiledEmbryos.AllHisRFP25CValidProfilesControlTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.HisRFP25CTimeSmoothedAvgAPProfiles.Control(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end

CompiledEmbryos.HisRFP25CTimeSmoothedAvgNarrowAPProfiles.Control = NaN(length(x), NumNarrowAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:NumNarrowAPbins
        BinData = NarrowDorsalProfiles(CompiledEmbryos.AllHisRFP25CValidProfilesControlTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.HisRFP25CTimeSmoothedAvgNarrowAPProfiles.Control(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end


%% Add TsetCsquishedTimes

CompiledEmbryos.AllTsetValidProfilesTestTF = CompiledEmbryos.IsNC14 &...
        CompiledEmbryos.TestSetEmbryos & ~isnan([CompiledEmbryos.TsetEmbryoTimes(:)].');% &...
if sum(CompiledEmbryos.AllTsetValidProfilesTestTF) > 0
y = CompiledEmbryos.TsetEmbryoTimes(CompiledEmbryos.AllTsetValidProfilesTestTF );
% ThreeMinDeltas = mink(y, 3);
% minx = floor(ThreeMinDeltas(end)/.1)*.1;
MaxT = ceil(max([CompiledEmbryos.TsetEmbryoTimes(:)]));
x = 0:1:MaxT;
sigma = sigma_time;
DiffMat = x.'-y;
[nn_idx, nn_D] = knnsearch(y.', x.', 'K', 25);
GaussianWeights = GetGaussianWeightMat(x, y, sigma, (width_sig_multiplier)*sigma);
keepx = ones(1, length(x), 'logical');
counts_within_1sigma = zeros(1, length(x));
counts_within_2sigma = zeros(1, length(x));
counts_within_3sigma = zeros(1, length(x));
counts_within_4sigma = zeros(1, length(x));
for x_index = 1:length(x)
    Num2SigmaPoints = sum(nn_D(x_index,:) <= 2*sigma_deltafc);
    
    if Num2SigmaPoints < min_2sigma_points
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    if sum(DiffMat(x_index, nn_idx(x_index,:)) > 0) < min_points_singleside | sum(DiffMat(x_index, nn_idx(x_index,:)) < 0) < min_points_singleside 
        keepx(x_index) = false;
        GaussianWeights(x_index,:) = NaN;
        continue
    end
    
    counts_within_1sigma(x_index) = sum(nn_D(x_index,:) <= sigma_deltafc);
    counts_within_2sigma(x_index) = sum(nn_D(x_index,:) <= 2*sigma_deltafc & nn_D(x_index,:) > sigma_deltafc);
    counts_within_3sigma(x_index) = sum(nn_D(x_index,:) <= 3*sigma_deltafc & nn_D(x_index,:) > 2*sigma_deltafc);
    counts_within_4sigma(x_index) = sum(nn_D(x_index,:) <= 4*sigma_deltafc & nn_D(x_index,:) > 3*sigma_deltafc);
    
    GaussianWeights(x_index, ~ismember(1:length(y), nn_idx(x_index,:))) = 0;
   
end
CompiledEmbryos.TsetSmoothedProfiles = {};
CompiledEmbryos.TsetSmoothedProfiles.TestCounts = {};
CompiledEmbryos.TsetSmoothedProfiles.TestCounts.counts_within_1sigma = counts_within_1sigma;
CompiledEmbryos.TsetSmoothedProfiles.TestCounts.counts_within_2sigma = counts_within_2sigma;
CompiledEmbryos.TsetSmoothedProfiles.TestCounts.counts_within_3sigma = counts_within_3sigma;
CompiledEmbryos.TsetSmoothedProfiles.TestCounts.counts_within_4sigma = counts_within_4sigma;


CompiledEmbryos.TsetSmoothedTimes= x;
CompiledEmbryos.TsetTimeSmoothedAvgAPProfiles = {};
CompiledEmbryos.TsetTimeSmoothedAvgAPProfiles.Test = NaN(length(x), NumAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:length(APbins)
        BinData = DorsalProfiles(CompiledEmbryos.AllTsetValidProfilesTestTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.TsetTimeSmoothedAvgAPProfiles.Test(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end

CompiledEmbryos.TsetTimeSmoothedAvgNarrowAPProfiles = {};
CompiledEmbryos.TsetTimeSmoothedAvgNarrowAPProfiles.Test = NaN(length(x), NumNarrowAPbins,  NChannels);

for ch_index = 2:NChannels % 2:5
    for RefBin = 1:NumNarrowAPbins
        BinData = NarrowDorsalProfiles(CompiledEmbryos.AllTsetValidProfilesTestTF,RefBin,ch_index);
        BinIsNaN = isnan(BinData);
        BinGaussianWeights = GaussianWeights(:,~BinIsNaN);
        WeightSums = sum(BinGaussianWeights, 2);
        BinData = BinData(~BinIsNaN);
        SmoothedBinProfile = ((BinGaussianWeights*BinData)./WeightSums).';
        CompiledEmbryos.TsetTimeSmoothedAvgNarrowAPProfiles.Test(:,RefBin,ch_index) = SmoothedBinProfile.';
    end
end
end

